#include <limits>
#include <cmath>
#include <vector>
#include <numeric>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include "Spectra/SymEigsSolver.h"
#include <nloptrAPI.h>
// [[Rcpp::depends(nloptr)]]
#include "utils.h"
#include "kernel.h"
#include "kriging.h"

double _nlopt_loss(unsigned n, const double* x, double* grad, void* data) {
    Kriging* kriging = static_cast<Kriging*>(data);
    std::size_t lengthscale_dim = kriging->get_lengthscale_dimension();
    Eigen::Map<Eigen::VectorXd> log_lengthscale(const_cast<double*>(x), lengthscale_dim);
    double log_nugget = kriging->is_interpolation() ? kriging->get_log_nugget()
                                                     : *(x + lengthscale_dim);
    return kriging->get_loss(log_lengthscale, log_nugget);
}

Kriging::Kriging(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, Kernel& Ker,
                 const bool& interpolation, const std::string& loss_type)
    : _n(X.rows()), _p(X.cols()), _X(X), _y(y), _Ker(Ker),
      _interpolation(interpolation), _L(_n) {
    _a.resize(_n);
    _b.resize(_n);
    _log_nugget = _interpolation ? std::log(1e-12) : std::log(1e-3);
    _set_default_lengthscale_bound();
    _validate_loss_type(loss_type);
    _set_y_tss();
}

void Kriging::_set_y_tss() {
    _y_tss = _y.squaredNorm() - std::pow(_y.sum(), 2) / _n;
}

void Kriging::_validate_loss_type(const std::string& loss_type) {
    if (loss_type == "NLLH") {
        _loss_type = loss_type;
    } else if (loss_type == "LOOCV") {
        Rcpp::stop("LOOCV loss is not implemented in this release.");
    } else {
        Rcpp::stop("Unknown loss type: \"%s\". Use \"NLLH\".", loss_type);
    }
}

bool Kriging::is_interpolation() const {
    return _interpolation;
}

bool Kriging::is_fitted() const {
    return _is_fitted;
}

void Kriging::_raise_error_if_not_fitted() const {
    if (!_is_fitted) Rcpp::stop("Model not fitted. Call fit() first.");
}

void Kriging::_invalidate_fit() {
    _loss = std::numeric_limits<double>::infinity();
    _is_fitted = false;
    _Kinv_pred.resize(0, 0);
}

void Kriging::_validate_nugget_bound() {
    if (_log_nugget_lower_bound > _log_nugget) {
        Rcpp::warning("Nugget lower bound (%g) exceeds current nugget (%g); clipping to current nugget.",
                      std::exp(_log_nugget_lower_bound), std::exp(_log_nugget));
        _log_nugget_lower_bound = _log_nugget;
    }
    if (_log_nugget_upper_bound < _log_nugget) {
        Rcpp::warning("Nugget upper bound (%g) is below current nugget (%g); clipping to current nugget.",
                      std::exp(_log_nugget_upper_bound), std::exp(_log_nugget));
        _log_nugget_upper_bound = _log_nugget;
    }
}

void Kriging::_validate_lengthscale_bound() {
    Eigen::VectorXd log_lengthscale = _Ker.get_log_lengthscale();
    std::size_t dim = log_lengthscale.size();
    for (std::size_t i = 0; i < dim; i++) {
        if (_log_lengthscale_lower_bound(i) > log_lengthscale(i)) {
            Rcpp::warning("Lengthscale lower bound (%g) exceeds current lengthscale (dim %d, %g); clipping to current lengthscale.",
                          std::exp(_log_lengthscale_lower_bound(i)), i + 1, std::exp(log_lengthscale(i)));
            _log_lengthscale_lower_bound(i) = log_lengthscale(i);
        }
        if (_log_lengthscale_upper_bound(i) < log_lengthscale(i)) {
            Rcpp::warning("Lengthscale upper bound (%g) is below current lengthscale (dim %d, %g); clipping to current lengthscale.",
                          std::exp(_log_lengthscale_upper_bound(i)), i + 1, std::exp(log_lengthscale(i)));
            _log_lengthscale_upper_bound(i) = log_lengthscale(i);
        }
    }
}

void Kriging::add_data(const Eigen::MatrixXd& Xn, const Eigen::VectorXd& yn) {
    _add_data(Xn, yn);
    _set_y_tss();
    _invalidate_fit();
}

void Kriging::_add_data(const Eigen::MatrixXd& Xn, const Eigen::VectorXd& yn) {
    // concatenate X and y with existing data
    std::size_t nn = Xn.rows();
    _X.conservativeResize(_n + nn, _p);
    _X.bottomRows(nn) = Xn;
    _y.conservativeResize(_n + nn);
    _y.tail(nn) = yn;
    // reset size for kriging parameters
    _n += nn;
    _a.resize(_n);
    _b.resize(_n);
    _L = Eigen::LLT<Eigen::MatrixXd>(_n);
    _Kinv_pred.resize(0, 0);
}


Rcpp::List Kriging::get_data() const {
    return Rcpp::List::create(Rcpp::Named("X") = _X, Rcpp::Named("y") = _y);
}

Eigen::MatrixXd Kriging::get_X() const {
    return _X;
}

Eigen::VectorXd Kriging::get_y() const {
    return _y;
}

std::size_t Kriging::get_datasize() const {
    return _n;
}

std::size_t Kriging::get_dimension() const {
    return _p;
}

std::size_t Kriging::get_lengthscale_dimension() const {
    return _Ker.get_log_lengthscale_dimension();
}

double Kriging::get_nllh() const {
    _raise_error_if_not_fitted();
    return _get_nllh();
}

double Kriging::get_loocv() const {
    Rcpp::stop("LOOCV loss is not implemented in this release.");
    return std::numeric_limits<double>::infinity();
}

double Kriging::get_loss() const {
    _raise_error_if_not_fitted();
    return _loss;
}

void Kriging::_set_default_lengthscale_bound() {
    Eigen::VectorXd log_lengthscale = _Ker.get_log_lengthscale();
    _log_lengthscale_lower_bound = (log_lengthscale.array() - std::log(1e3)).matrix();
    _log_lengthscale_upper_bound = (log_lengthscale.array() + std::log(1e3)).matrix();
}

void Kriging::_optimize_hyperparameters() {
    std::size_t lengthscale_dim = get_lengthscale_dimension();
    Eigen::VectorXd x(lengthscale_dim + 1), lb(lengthscale_dim + 1), ub(lengthscale_dim + 1);
    x.head(lengthscale_dim) = _Ker.get_log_lengthscale();
    x(lengthscale_dim) = _log_nugget;
    lb.head(lengthscale_dim) = _log_lengthscale_lower_bound;
    lb(lengthscale_dim) = _log_nugget_lower_bound;
    ub.head(lengthscale_dim) = _log_lengthscale_upper_bound;
    ub(lengthscale_dim) = _log_nugget_upper_bound;

    std::size_t opt_dim = _interpolation ? lengthscale_dim : lengthscale_dim + 1;
    nlopt_opt optimizer = _nlopt_init(_nlopt_algorithm, opt_dim);
    if (_nlopt_local_algorithm != "") {
        nlopt_opt local_optimizer = _nlopt_init(_nlopt_local_algorithm, opt_dim);
        nlopt_set_xtol_rel(local_optimizer, 1e-4);
        nlopt_set_local_optimizer(optimizer, local_optimizer);
        nlopt_set_population(optimizer, _nlopt_population);
        nlopt_destroy(local_optimizer);
    }
    nlopt_set_lower_bounds(optimizer, lb.data());
    nlopt_set_upper_bounds(optimizer, ub.data());
    nlopt_set_min_objective(optimizer, _nlopt_loss, this);
    nlopt_set_maxeval(optimizer, _nlopt_maxeval);

    double loss_min;
    nlopt_optimize(optimizer, x.data(), &loss_min);
    nlopt_destroy(optimizer);
}


void Kriging::set_lengthscale(const Eigen::VectorXd& lengthscale) {
    _Ker.set_lengthscale(lengthscale);
    _validate_lengthscale_bound();
    _invalidate_fit();
}

Eigen::VectorXd Kriging::get_lengthscale() const {
    return _Ker.get_lengthscale();
}

void Kriging::relax_lengthscale_constraint() {
    _Ker.relax_lengthscale_constraint();
    _set_default_lengthscale_bound();
}

void Kriging::set_nugget(const double& nugget) {
    if (nugget < 0)
        Rcpp::stop("Nugget must be non-negative.");
    _log_nugget = std::log(std::max(nugget, STABILITY_THRESHOLD));
    _validate_nugget_bound();
    _invalidate_fit();
}

void Kriging::set_nugget_bound(const double& nugget_lower_bound,
                               const double& nugget_upper_bound) {
    if (nugget_lower_bound < 0 || nugget_upper_bound < 0)
        Rcpp::stop("Nugget bounds must be non-negative.");
    if (nugget_lower_bound > nugget_upper_bound)
        Rcpp::stop("Nugget lower bound must be less than or equal to upper bound.");
    _log_nugget_lower_bound = std::log(std::max(nugget_lower_bound, STABILITY_THRESHOLD));
    _log_nugget_upper_bound = std::log(std::max(nugget_upper_bound, STABILITY_THRESHOLD));
    _validate_nugget_bound();
}

double Kriging::get_log_nugget() const {
    return _log_nugget;
}

void Kriging::set_epsilon(const double& epsilon) {
    _epsilon = epsilon;
}

void Kriging::set_nlopt_algorithm(const std::string& nlopt_algorithm) {
    _nlopt_algorithm = nlopt_algorithm;
}

void Kriging::set_nlopt_local_algorithm(const std::string& nlopt_local_algorithm) {
    _nlopt_local_algorithm = nlopt_local_algorithm;
}

void Kriging::set_nlopt_maxeval(const std::size_t& nlopt_maxeval) {
    _nlopt_maxeval = nlopt_maxeval;
}

void Kriging::set_nlopt_population(const std::size_t& nlopt_population) {
    _nlopt_population = nlopt_population;
}

std::string Kriging::get_loss_type() const {
    return _loss_type;
}

void Kriging::set_lengthscale_bound(const Eigen::VectorXd& lengthscale_lower_bound,
                                    const Eigen::VectorXd& lengthscale_upper_bound) {
    if (lengthscale_lower_bound.size() != get_lengthscale_dimension() || lengthscale_upper_bound.size() != get_lengthscale_dimension())
        Rcpp::stop("Lengthscale bound dimension does not match kriging dimension.");
    if ((lengthscale_lower_bound.array() <= 0).any() || (lengthscale_upper_bound.array() <= 0).any())
        Rcpp::stop("Lengthscale bounds must be positive.");
    if ((lengthscale_lower_bound.array() > lengthscale_upper_bound.array()).any())
        Rcpp::stop("Lengthscale lower bound must be less than or equal to upper bound.");
    _log_lengthscale_lower_bound = lengthscale_lower_bound.array().log().matrix();
    _log_lengthscale_upper_bound = lengthscale_upper_bound.array().log().matrix();
    _validate_lengthscale_bound();
}

void Kriging::fit(bool optimize) {
    if (optimize) _optimize_hyperparameters();
    _set_hyperparameters();
}

double Kriging::get_loss(const Eigen::VectorXd& log_lengthscale, const double& log_nugget) {
    if (!_set_hyperparameters(log_lengthscale, log_nugget)) {
        return std::numeric_limits<double>::infinity();
    }
    double loss = (_loss_type == "NLLH") ? _get_nllh() : _get_loocv();
    // Track best: update kernel/nugget only on improvement
    if (loss < _loss) {
        _Ker.set_log_lengthscale(log_lengthscale);
        _log_nugget = log_nugget;
        _loss = loss;
    }
    return loss;
}

void Kriging::_set_hyperparameters() {
    if (!_set_hyperparameters(_Ker.get_log_lengthscale(), _log_nugget)) {
        Rcpp::stop(
            "Cholesky decomposition failed. "
            "The current lengthscale and nugget yield a non-positive definite kernel matrix. "
            "Try using smaller lengthscale or larger nugget."
        );
    }
    if (_interpolation) _update_prediction_inverse_neumann();
    _loss = (_loss_type == "NLLH") ? _get_nllh() : _get_loocv();
    _is_fitted = true;
}

double Kriging::get_mu() const {
    _raise_error_if_not_fitted();
    return _mu;
}

double Kriging::get_nu2() const {
    _raise_error_if_not_fitted();
    return _nu2;
}

double Kriging::get_sigma2() const {
    _raise_error_if_not_fitted();
    return _interpolation ? 0 : std::exp(_log_nugget) * _nu2;
}

double Kriging::get_nugget() const {
    return std::exp(_log_nugget);
}

void Kriging::_update_prediction_inverse_neumann(const std::size_t& max_iter) {
    Eigen::MatrixXd A = _L.solve(Eigen::MatrixXd::Identity(_n, _n));
    _Kinv_pred = A;
    double nugget = std::exp(_log_nugget);
    if (nugget > 0.0) {
        Eigen::MatrixXd term = A;
        for (std::size_t iter = 0; iter < max_iter; ++iter) {
            term = nugget * (A * term);
            _Kinv_pred += term;
            if (term.norm() <= STABILITY_THRESHOLD * std::max(1.0, _Kinv_pred.norm())) break;
        }
    }
    _Kinv_pred = 0.5 * (_Kinv_pred + _Kinv_pred.transpose()).eval();
}

Eigen::VectorXd Kriging::predict(const Eigen::VectorXd& xnew) const {
    _raise_error_if_not_fitted();
    double mean, sd;
    _predict(xnew, mean, sd);
    return Eigen::Vector2d(mean, sd);
}

Rcpp::List Kriging::predict(const Eigen::MatrixXd& Xnew) const {
    _raise_error_if_not_fitted();
    std::size_t m(Xnew.rows());
    Eigen::VectorXd vec_mean(m), vec_sd(m);
    double mean, sd;
    for (std::size_t i = 0; i < m; i++) {
        _predict(Xnew.row(i), mean, sd);
        vec_mean(i) = mean;
        vec_sd(i) = sd;
    }
    return Rcpp::List::create(Rcpp::Named("mean") = vec_mean, Rcpp::Named("sd") = vec_sd);
}

Rcpp::List Kriging::predict_neumann(const Eigen::MatrixXd& Xnew) const {
    _raise_error_if_not_fitted();
    if (!_interpolation) return predict(Xnew);
    return _predict_neumann(Xnew);
}

// =============================================================================
// OrdinaryKriging
// =============================================================================

bool OrdinaryKriging::_set_hyperparameters(const Eigen::VectorXd& log_lengthscale, const double& log_nugget) {
    Eigen::MatrixXd R = _Ker.compute(_X, log_lengthscale);
    R.diagonal().array() += std::exp(log_nugget);
    _L = R.llt();
    if (_L.info() != 0) return false;
    _a = _L.matrixL().solve(Eigen::VectorXd::Ones(_n));
    _b = _L.matrixL().solve(_y);
    _mu = _a.dot(_b) / _a.dot(_a);
    _nu2 = std::max(1.0 / (_n - 1) * (_b.dot(_b) - _mu * _mu * _a.dot(_a)), STABILITY_THRESHOLD);
    return true;
}

double OrdinaryKriging::_get_nllh() const {
    double nllh = (_n - 1) * std::log(_nu2) 
                + 2 * _L.matrixLLT().diagonal().array().log().sum() 
                + std::log(_a.dot(_a));
    if (_interpolation) {
        Eigen::VectorXd ym = _b - _mu * _a;
        Eigen::VectorXd pred = (_mu + (_L.matrixL() * ym - _L.matrixU().solve(ym) * std::exp(_log_nugget)).array()).matrix();
        double penalty = 2 * (_y - pred).squaredNorm() / (_epsilon * _y_tss);
        nllh = nllh / _n + penalty;
    }
    return nllh;
}

double OrdinaryKriging::_get_loocv() const {
    // NOT IMPLEMENTED YET
    return std::numeric_limits<double>::infinity();
}

void OrdinaryKriging::_predict(const Eigen::VectorXd& xnew, double& mean, double& sd) const {
    Eigen::VectorXd r = _Ker.compute(xnew, _X);
    Eigen::VectorXd d = _L.matrixL().solve(r);
    mean = _mu + d.dot(_b - _mu * _a);
    sd = std::sqrt(_nu2) * std::sqrt(
        std::max(1 - d.dot(d) + std::pow(1 - d.dot(_a), 2) / _a.dot(_a), 0.0)
    );
}

Rcpp::List OrdinaryKriging::_predict_neumann(const Eigen::MatrixXd& Xnew) const {
    std::size_t m = Xnew.rows();
    Eigen::VectorXd vec_mean(m), vec_sd(m);
    Eigen::VectorXd pred_inv_one = _Kinv_pred * Eigen::VectorXd::Ones(_n);
    double pred_one_inv_one = pred_inv_one.sum();
    Eigen::VectorXd pred_inv_y = _Kinv_pred * _y;
    for (std::size_t i = 0; i < m; ++i) {
        Eigen::VectorXd xnew = Xnew.row(i).transpose();
        Eigen::VectorXd r = _Ker.compute(xnew, _X);
        Eigen::VectorXd d = _Kinv_pred * r;
        double one_term = r.dot(pred_inv_one);
        vec_mean(i) = _mu + r.dot(pred_inv_y - _mu * pred_inv_one);
        vec_sd(i) = std::sqrt(_nu2) * std::sqrt(
            std::max(1.0 - r.dot(d) + std::pow(1.0 - one_term, 2) / pred_one_inv_one, 0.0)
        );
    }
    return Rcpp::List::create(Rcpp::Named("mean") = vec_mean, Rcpp::Named("sd") = vec_sd);
}

// =============================================================================
// LimitKriging
// =============================================================================

double LimitKriging::_get_loocv() const {
    // NOT IMPLEMENTED YET
    return std::numeric_limits<double>::infinity();
}

void LimitKriging::_predict(const Eigen::VectorXd& xnew, double& mean, double& sd) const {
    Eigen::VectorXd r = _Ker.compute(xnew, _X);
    Eigen::VectorXd d = _L.matrixL().solve(r);
    mean = d.dot(_b) / d.dot(_a);
    sd = std::sqrt(_nu2) * std::sqrt(
        std::max(1 - d.dot(d) + d.dot(d) * std::pow(1 - d.dot(_a), 2) / std::pow(d.dot(_a), 2), 0.0)
    );
}

Rcpp::List LimitKriging::_predict_neumann(const Eigen::MatrixXd& Xnew) const {
    std::size_t m = Xnew.rows();
    Eigen::VectorXd vec_mean(m), vec_sd(m);
    Eigen::VectorXd pred_inv_one = _Kinv_pred * Eigen::VectorXd::Ones(_n);
    Eigen::VectorXd pred_inv_y = _Kinv_pred * _y;
    for (std::size_t i = 0; i < m; ++i) {
        Eigen::VectorXd xnew = Xnew.row(i).transpose();
        Eigen::VectorXd r = _Ker.compute(xnew, _X);
        Eigen::VectorXd d = _Kinv_pred * r;
        double one_term = r.dot(pred_inv_one);
        double quad_term = r.dot(d);
        vec_mean(i) = r.dot(pred_inv_y) / one_term;
        vec_sd(i) = std::sqrt(_nu2) * std::sqrt(
            std::max(1.0 - quad_term + quad_term * std::pow(1.0 - one_term, 2) / std::pow(one_term, 2), 0.0)
        );
    }
    return Rcpp::List::create(Rcpp::Named("mean") = vec_mean, Rcpp::Named("sd") = vec_sd);
}

// =============================================================================
// RationalKriging
// =============================================================================

RationalKriging::RationalKriging(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, 
                                 Kernel& Ker, const bool& interpolation, 
                                 const std::string& loss_type)
    : Kriging(X, y, Ker, interpolation, loss_type) {
    _c.resize(_n);
    _s.resize(_n);
}

void RationalKriging::_add_data(const Eigen::MatrixXd& Xn, const Eigen::VectorXd& yn) {
    Kriging::_add_data(Xn, yn);
    _c.resize(_n);
    _s.resize(_n);
}

bool RationalKriging::_set_hyperparameters(const Eigen::VectorXd& log_lengthscale, const double& log_nugget) {
    Eigen::MatrixXd R = _Ker.compute(_X, log_lengthscale);
    R.diagonal().array() += std::exp(log_nugget);
    // get eigen vectors / values (rank from smallest to largest)
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigSolver(R);
    // compute gamma and c via binary search
    Eigen::VectorXd eigOne = eigSolver.eigenvectors().transpose() * Eigen::VectorXd::Ones(_n);
    double delta = eigSolver.eigenvalues().maxCoeff() / _n;
    double gamma, gamma_lb = 0.0, gamma_ub = 1.0, eps = 1e-4;
    _c = eigSolver.eigenvectors() * (
        eigOne.array() / ((1 - gamma_lb) * eigSolver.eigenvalues().array() + gamma_lb)
    ).matrix();
    if (_c.minCoeff() > delta) {
        gamma = gamma_lb;
    } else {
        while (gamma_ub - gamma_lb > eps) {
            gamma = 0.5 * (gamma_lb + gamma_ub);
            _c = eigSolver.eigenvectors() * (
                eigOne.array() / ((1 - gamma) * eigSolver.eigenvalues().array() + gamma)
            ).matrix();
            if (_c.minCoeff() < delta) {
                gamma_lb = gamma;
            } else {
                gamma_ub = gamma;
            }
        }
        gamma = gamma_ub;
        _c = eigSolver.eigenvectors() * (
            eigOne.array() / ((1 - gamma) * eigSolver.eigenvalues().array() + gamma)
        ).matrix();
    }
    _c = _c.cwiseMax(0).normalized();
    // compute parameters
    _s = R * _c;
    _L = R.llt();
    if (_L.info() != 0) return false;
    _a = _L.matrixL().solve(_s);
    _b = _L.matrixL().solve(_s.cwiseProduct(_y));
    _mu = _a.dot(_b) / _a.dot(_a);
    _nu2 = std::max(1.0 / (_n - 1) * (_b - _mu * _a).squaredNorm(), STABILITY_THRESHOLD);
    return true;
}

double RationalKriging::_get_nllh() const {
    double nllh = (_n - 1) * std::log(_nu2) 
                + 2 * _L.matrixLLT().diagonal().array().log().sum() 
                - 2 * _s.array().log().sum() 
                + std::log(_c.dot(_s));
    if (_interpolation) {
        Eigen::VectorXd pred = (_L.matrixL() * _b - _L.matrixU().solve(_b) * std::exp(_log_nugget)).cwiseQuotient(_s);
        double penalty = 2 * (_y - pred).squaredNorm() / (_epsilon * _y_tss);
        nllh = nllh / _n + penalty;
    }
    return nllh;
}

double RationalKriging::_get_loocv() const {
    // NOT IMPLEMENTED YET
    return std::numeric_limits<double>::infinity();
}

Eigen::VectorXd RationalKriging::get_c() const {
    _raise_error_if_not_fitted();
    return _c;
}

void RationalKriging::_predict(const Eigen::VectorXd& xnew, double& mean, double& sd) const {
    Eigen::VectorXd r = _Ker.compute(xnew, _X);
    Eigen::VectorXd d = _L.matrixL().solve(r);
    mean = d.dot(_b) / d.dot(_a);
    sd = std::sqrt(_nu2) * std::sqrt(std::max(1 - d.dot(d), 0.0)) / d.dot(_a);
}

Rcpp::List RationalKriging::_predict_neumann(const Eigen::MatrixXd& Xnew) const {
    std::size_t m = Xnew.rows();
    Eigen::VectorXd vec_mean(m), vec_sd(m);
    Eigen::VectorXd s_pred = _s - std::exp(_log_nugget) * _c;
    Eigen::VectorXd pred_inv_sy = _Kinv_pred * s_pred.cwiseProduct(_y);
    for (std::size_t i = 0; i < m; ++i) {
        Eigen::VectorXd xnew = Xnew.row(i).transpose();
        Eigen::VectorXd r = _Ker.compute(xnew, _X);
        Eigen::VectorXd d = _Kinv_pred * r;
        double sx = r.dot(_c);
        vec_mean(i) = r.dot(pred_inv_sy) / sx;
        vec_sd(i) = std::sqrt(_nu2) * std::sqrt(std::max(1.0 - r.dot(d), 0.0)) / sx;
    }
    return Rcpp::List::create(Rcpp::Named("mean") = vec_mean, Rcpp::Named("sd") = vec_sd);
}

// =============================================================================
// GeneralizedRationalKriging
// =============================================================================

GeneralizedRationalKriging::GeneralizedRationalKriging(const Eigen::MatrixXd& X, 
                                                       const Eigen::VectorXd& y, 
                                                       Kernel& Ker, const bool& interpolation, 
                                                       const std::string& loss_type)
    : Kriging(X, y, Ker, interpolation, loss_type) {
    _c.resize(_n + 1);
    _s.resize(_n);
}

void GeneralizedRationalKriging::_add_data(const Eigen::MatrixXd& Xn, const Eigen::VectorXd& yn) {
    Kriging::_add_data(Xn, yn);
    _c.resize(_n + 1);
    _s.resize(_n);
}

bool GeneralizedRationalKriging::_set_hyperparameters(const Eigen::VectorXd& log_lengthscale, const double& log_nugget) {
    Eigen::MatrixXd R = _Ker.compute(_X, log_lengthscale);
    R.diagonal().array() += std::exp(log_nugget);
    _L = R.llt();
    if (_L.info() != 0) return false;
    // compute (c0, c')
    Eigen::MatrixXd M(_n + 1, _n + 1);
    M(0, 0) = _L.matrixL().solve(Eigen::VectorXd::Ones(_n)).squaredNorm();
    M.bottomLeftCorner(_n, 1) = Eigen::VectorXd::Ones(_n);
    M.topRightCorner(1, _n) = Eigen::VectorXd::Ones(_n);
    M.bottomRightCorner(_n, _n) = R;
    Spectra::DenseSymMatProd<double> op(M);
    Spectra::SymEigsSolver<Spectra::DenseSymMatProd<double>> eigs(op, 1, (_n + 1 < 20) ? _n + 1 : 20);
    eigs.init();
    eigs.compute(Spectra::SortRule::LargestAlge);
    _c = eigs.eigenvectors().col(0).cwiseAbs();
    _s = (_c(0) + (R * _c.tail(_n)).array()).matrix();
    _a = _L.matrixL().solve(_s);
    _b = _L.matrixL().solve(_s.cwiseProduct(_y));
    _mu = _a.dot(_b) / _a.dot(_a);
    _nu2 = std::max(1.0 / (_n - 1) * (_b - _mu * _a).squaredNorm(), STABILITY_THRESHOLD);
    return true;
}

double GeneralizedRationalKriging::_get_nllh() const {
    double nllh = (_n - 1) * std::log(_nu2) 
                + 2 * _L.matrixLLT().diagonal().array().log().sum() 
                - 2 * _s.array().log().sum() 
                + std::log(_a.dot(_a));
    if (_interpolation) {
        Eigen::VectorXd ym = _b - _mu * _a;
        Eigen::VectorXd pred = (
            _mu + (_L.matrixL() * ym - _L.matrixU().solve(ym) * std::exp(_log_nugget)).cwiseQuotient(_s).array()
        ).matrix();
        double penalty = 2 * (_y - pred).squaredNorm() / (_epsilon * _y_tss);
        nllh = nllh / _n + penalty;
    }
    return nllh;
}

double GeneralizedRationalKriging::_get_loocv() const {
    // NOT IMPLEMENTED YET
    return std::numeric_limits<double>::infinity();
}

double GeneralizedRationalKriging::get_c0() const {
    _raise_error_if_not_fitted();
    return _c(0);
}

Eigen::VectorXd GeneralizedRationalKriging::get_c() const {
    _raise_error_if_not_fitted();
    return _c.tail(_n);
}

void GeneralizedRationalKriging::_predict(const Eigen::VectorXd& xnew, double& mean, double& sd) const {
    Eigen::VectorXd r = _Ker.compute(xnew, _X);
    Eigen::VectorXd d = _L.matrixL().solve(r);
    double sx = _c(0) + r.dot(_c.tail(_n));
    double tx = sx - d.dot(_a);
    mean = (_mu * tx + d.dot(_b)) / sx;
    sd = std::sqrt(_nu2) * std::sqrt(
        std::max(1 - d.dot(d) + std::pow(tx, 2) / _a.dot(_a), 0.0)
    ) / sx;
}

Rcpp::List GeneralizedRationalKriging::_predict_neumann(const Eigen::MatrixXd& Xnew) const {
    std::size_t m = Xnew.rows();
    Eigen::VectorXd vec_mean(m), vec_sd(m);
    Eigen::VectorXd pred_inv_one = _Kinv_pred * Eigen::VectorXd::Ones(_n);
    double pred_one_inv_one = pred_inv_one.sum();
    Eigen::VectorXd s_pred = _s - std::exp(_log_nugget) * _c.tail(_n);
    Eigen::VectorXd pred_inv_sy = _Kinv_pred * s_pred.cwiseProduct(_y);
    for (std::size_t i = 0; i < m; ++i) {
        Eigen::VectorXd xnew = Xnew.row(i).transpose();
        Eigen::VectorXd r = _Ker.compute(xnew, _X);
        Eigen::VectorXd d = _Kinv_pred * r;
        double sx = _c(0) + r.dot(_c.tail(_n));
        double tx = _c(0) * (1.0 - r.dot(pred_inv_one));
        vec_mean(i) = (_mu * tx + r.dot(pred_inv_sy)) / sx;
        vec_sd(i) = std::sqrt(_nu2) * std::sqrt(
            std::max(1.0 - r.dot(d) + std::pow(tx, 2) / pred_one_inv_one, 0.0)
        ) / sx;
    }
    return Rcpp::List::create(Rcpp::Named("mean") = vec_mean, Rcpp::Named("sd") = vec_sd);
}

// =============================================================================
// UniversalKriging
// =============================================================================

UniversalKriging::UniversalKriging(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, 
                                   Kernel& Ker, const bool& interpolation, 
                                   const std::size_t& pb, Rcpp::Function bfunc, 
                                   const std::string& loss_type)
    : Kriging(X, y, Ker, interpolation, loss_type), _pb(pb), _bfunc(bfunc), _LQ(_n, _pb) {
    _beta.resize(_pb);
    _B.resize(_n, _pb);
    for (std::size_t i = 0; i < _n; i++) {
        std::vector<double> basis = Rcpp::as<std::vector<double>>(_bfunc(_X.row(i)));
        _B.row(i) = Eigen::Map<Eigen::VectorXd>(&basis[0], _pb);
    }
}

void UniversalKriging::_add_data(const Eigen::MatrixXd& Xn, const Eigen::VectorXd& yn) {
    std::size_t nn = Xn.rows();
    Eigen::MatrixXd Bn(nn, _pb);
    for (std::size_t i = 0; i < nn; i++) {
        std::vector<double> basis = Rcpp::as<std::vector<double>>(_bfunc(Xn.row(i)));
        Bn.row(i) = Eigen::Map<Eigen::VectorXd>(&basis[0], _pb);
    }
    Eigen::Map<Eigen::MatrixXd> Bo(_B.data(), _n, _pb);
    _B.resize(_n + nn, _pb);
    _B.topRows(_n) = Bo;
    _B.bottomRows(nn) = Bn;
    Kriging::_add_data(Xn, yn);
    _LQ = Eigen::ColPivHouseholderQR<Eigen::MatrixXd>(_n, _pb);
}

bool UniversalKriging::_set_hyperparameters(const Eigen::VectorXd& log_lengthscale, const double& log_nugget) {
    Eigen::MatrixXd R = _Ker.compute(_X, log_lengthscale);
    R.diagonal().array() += std::exp(log_nugget);
    _L = R.llt();
    if (_L.info() != 0) return false;
    _LQ = Eigen::ColPivHouseholderQR<Eigen::MatrixXd>(_L.matrixL().solve(_B));
    _b = _L.matrixL().solve(_y);
    _beta = _LQ.solve(_b);
    Eigen::VectorXd mu = _B * _beta;
    _a = _L.matrixL().solve(mu);
    _nu2 = std::max(1.0 / (_n - _pb) * (_b - _a).squaredNorm(), STABILITY_THRESHOLD);
    return true;
}

double UniversalKriging::_get_nllh() const {
    double nllh = (_n - _pb) * std::log(_nu2) 
                + 2 * _L.matrixLLT().diagonal().array().log().sum() 
                + 2 * _LQ.logAbsDeterminant();
    if (_interpolation) {
        Eigen::VectorXd mu = _B * _beta;
        Eigen::VectorXd ym = _b - _a;
        Eigen::VectorXd pred = mu + (_L.matrixL() * ym - _L.matrixU().solve(ym) * std::exp(_log_nugget));
        double penalty = 2 * (_y - pred).squaredNorm() / (
            _epsilon * std::min(_y_tss, (_y - mu).squaredNorm() - std::pow((_y - mu).sum(), 2) / _n)
        );
        nllh = nllh / _n + penalty;
    }
    return nllh;
}

double UniversalKriging::_get_loocv() const {
    // NOT IMPLEMENTED YET
    return std::numeric_limits<double>::infinity();
}

Eigen::VectorXd UniversalKriging::get_beta() const {
    _raise_error_if_not_fitted();
    return _beta;
}

void UniversalKriging::_predict(const Eigen::VectorXd& xnew, double& mean, double& sd) const {
    std::vector<double> basis = Rcpp::as<std::vector<double>>(_bfunc(xnew));
    Eigen::Map<Eigen::VectorXd> bm(&basis[0], _pb);
    Eigen::VectorXd r = _Ker.compute(xnew, _X);
    Eigen::VectorXd d = _L.matrixL().solve(r);
    Eigen::VectorXd h = _LQ.colsPermutation().transpose() * bm 
                      - _LQ.matrixQR().triangularView<Eigen::Upper>().transpose() 
                      * (_LQ.householderQ().transpose() * d);
    mean = bm.dot(_beta) + d.dot(_b - _a);
    sd = std::sqrt(_nu2) * std::sqrt(
        std::max(1 - d.dot(d) 
               + (_LQ.matrixQR().topLeftCorner(_pb, _pb).triangularView<Eigen::Upper>()
                  .transpose().solve(h)).squaredNorm(), 0.0)
    );
}

Rcpp::List UniversalKriging::_predict_neumann(const Eigen::MatrixXd& Xnew) const {
    std::size_t m = Xnew.rows();
    Eigen::VectorXd vec_mean(m), vec_sd(m);
    Eigen::VectorXd mu = _B * _beta;
    Eigen::VectorXd pred_resid = _Kinv_pred * (_y - mu);
    Eigen::MatrixXd pred_inv_B = _Kinv_pred * _B;
    Eigen::MatrixXd pred_Bt_inv_B = _B.transpose() * pred_inv_B;
    pred_Bt_inv_B = 0.5 * (pred_Bt_inv_B + pred_Bt_inv_B.transpose()).eval();
    for (std::size_t i = 0; i < m; ++i) {
        Eigen::VectorXd xnew = Xnew.row(i).transpose();
        std::vector<double> basis = Rcpp::as<std::vector<double>>(_bfunc(xnew));
        Eigen::Map<Eigen::VectorXd> bm(&basis[0], _pb);
        Eigen::VectorXd r = _Ker.compute(xnew, _X);
        Eigen::VectorXd d = _Kinv_pred * r;
        Eigen::VectorXd h = bm - pred_inv_B.transpose() * r;
        vec_mean(i) = bm.dot(_beta) + r.dot(pred_resid);
        vec_sd(i) = std::sqrt(_nu2) * std::sqrt(
            std::max(1.0 - r.dot(d) + h.dot(pred_Bt_inv_B.ldlt().solve(h)), 0.0)
        );
    }
    return Rcpp::List::create(Rcpp::Named("mean") = vec_mean, Rcpp::Named("sd") = vec_sd);
}
