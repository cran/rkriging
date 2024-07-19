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

double nlopt_nllh(unsigned n, const double* x, double* grad, void* data) {
    // cast void* pointer to a pointer of type Kriging
    Kriging* kriging = (Kriging*) data;
    std::size_t lengthscale_dim = kriging->get_lengthscale_dimension();
    Eigen::Map<Eigen::VectorXd> log_lengthscale(const_cast<double*>(x), lengthscale_dim);
    double nugget = kriging->is_interpolation() ? 1e-6 : std::exp(*(x+lengthscale_dim)); 
    return kriging->get_nllh(log_lengthscale, nugget);
}

Kriging::Kriging(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, Kernel& Ker, const bool& interpolation): n_(X.rows()), p_(X.cols()), X_(X), y_(y), Ker_(Ker), interpolation_(interpolation), L_(n_) {
    a_.resize(n_);
    b_.resize(n_);
    // L_ = Eigen::LLT<Eigen::MatrixXd>(n_);
    y_tss_ = y_.squaredNorm() - std::pow(y_.sum(),2)/n_;
    nugget_ = interpolation_ ? 1e-6 : 1e-3;
    nllh_ = std::numeric_limits<double>::infinity();
}

bool Kriging::is_interpolation() {
    return interpolation_;
}

void Kriging::add_data(const Eigen::MatrixXd& Xn, const Eigen::VectorXd& yn) {
    std::size_t nn = Xn.rows();
    Eigen::Map<Eigen::MatrixXd> Xo(X_.data(), n_, p_);
    X_.resize(n_+nn, p_);
    X_.topRows(n_) = Xo;
    X_.bottomRows(nn) = Xn;
    Eigen::Map<Eigen::VectorXd> yo(y_.data(), n_);
    y_.resize(n_+nn);
    y_.head(n_) = yo;
    y_.tail(nn) = yn;
    n_ += nn;
    a_.resize(n_);
    b_.resize(n_);
    L_ = Eigen::LLT<Eigen::MatrixXd>(n_);
    y_tss_ = y_.squaredNorm() - std::pow(y_.sum(),2)/n_;
    nllh_ = std::numeric_limits<double>::infinity();
}

Rcpp::List Kriging::get_data() {
    return Rcpp::List::create(Rcpp::Named("X") = X_, Rcpp::Named("y") = y_);
}

Eigen::MatrixXd Kriging::get_X() {
    return X_;
}

Eigen::VectorXd Kriging::get_y() {
    return y_;
}

std::size_t Kriging::get_datasize() {
    return n_;
}

std::size_t Kriging::get_dimension() {
    return p_;
}

std::size_t Kriging::get_lengthscale_dimension() {
    return Ker_.get_log_lengthscale_dimension();
}

double Kriging::get_nllh() {
    return nllh_;
}

void Kriging::fit_hyperparameters(const Eigen::VectorXd& lengthscale_lower_bound, const Eigen::VectorXd& lengthscale_upper_bound) {
    std::size_t lengthscale_dim = get_lengthscale_dimension();
    Eigen::VectorXd x(lengthscale_dim+1), lb(lengthscale_dim+1), ub(lengthscale_dim+1);
    x.head(lengthscale_dim) = Ker_.get_log_lengthscale();
    x(lengthscale_dim) = std::log(nugget_);
    lb.head(lengthscale_dim) = lengthscale_lower_bound.array().log().matrix();
    lb(lengthscale_dim) = std::log(1e-6);
    ub.head(lengthscale_dim) = lengthscale_upper_bound.array().log().matrix();
    ub(lengthscale_dim) = 0;
    
    std::size_t opt_dim = interpolation_ ? lengthscale_dim : lengthscale_dim+1;
    nlopt_opt optimizer = nlopt_init(nlopt_algorithm_, opt_dim);
    if (nlopt_local_algorithm_ != "") {
        nlopt_opt local_optimizer = nlopt_init(nlopt_local_algorithm_, opt_dim);
        nlopt_set_xtol_rel(local_optimizer, 1e-4);
        nlopt_set_local_optimizer(optimizer, local_optimizer);
        nlopt_set_population(optimizer, (opt_dim+1 > 10) ? opt_dim+1 : 10);
        nlopt_destroy(local_optimizer);
    }
    nlopt_set_lower_bounds(optimizer, lb.data());
    nlopt_set_upper_bounds(optimizer, ub.data());
    nlopt_set_min_objective(optimizer, nlopt_nllh, this);
    nlopt_set_maxeval(optimizer, nlopt_maxeval_);

    double nllh_min;
    nlopt_optimize(optimizer, x.data(), &nllh_min);
    nlopt_destroy(optimizer);

    Ker_.set_log_lengthscale(x.head(lengthscale_dim));
    if (!interpolation_) nugget_ = std::exp(x(lengthscale_dim));
}

void Kriging::set_lengthscale(const Eigen::VectorXd& lengthscale) {
    Ker_.set_lengthscale(lengthscale);
    nllh_ = std::numeric_limits<double>::infinity();
}

Eigen::VectorXd Kriging::get_lengthscale() {
    return Ker_.get_lengthscale();
}

void Kriging::relax_lengthscale_constraint() {
    Ker_.relax_lengthscale_constraint();
}

void Kriging::set_nugget(const double& nugget) {
    nugget_ = nugget;
}

void Kriging::set_epsilon(const double& epsilon) {
    epsilon_ = epsilon;
}

void Kriging::set_nlopt_algorithm(const std::string& nlopt_algorithm) {
    nlopt_algorithm_ = nlopt_algorithm;
}

void Kriging::set_nlopt_local_algorithm(const std::string& nlopt_local_algorithm) {
    nlopt_local_algorithm_ = nlopt_local_algorithm;
}

void Kriging::set_nlopt_maxeval(const std::size_t& nlopt_maxeval) {
    nlopt_maxeval_ = nlopt_maxeval;
}

void Kriging::fit() {
    Eigen::VectorXd log_lengthscale = Ker_.get_log_lengthscale();
    fit((log_lengthscale.array()-std::log(1e3)).exp().matrix(), (log_lengthscale.array()+std::log(1e3)).exp().matrix());
}

void Kriging::fit(const Eigen::VectorXd& lengthscale_lower_bound, const Eigen::VectorXd& lengthscale_upper_bound) {
    fit_hyperparameters(lengthscale_lower_bound, lengthscale_upper_bound);
    set_kriging_parameters();
}

double Kriging::get_mu() {
    return mu_;
}

double Kriging::get_nu2() {
    return nu2_;
}

double Kriging::get_sigma2() {
    return interpolation_ ? 0 : nugget_ * nu2_;
}

Eigen::VectorXd Kriging::predict(const Eigen::VectorXd& xnew) {
    double mean, sd;
    predict(xnew, mean, sd);
    return Eigen::Vector2d(mean, sd);
}

Rcpp::List Kriging::predict(const Eigen::MatrixXd& Xnew) {
    std::size_t m(Xnew.rows());
    Eigen::VectorXd vec_mean(m), vec_sd(m);
    double mean, sd;
    for (std::size_t i = 0; i < m; i++) {
        predict(Xnew.row(i), mean, sd);
        vec_mean(i) = mean;
        vec_sd(i) = sd;
    }
    return Rcpp::List::create(Rcpp::Named("mean") = vec_mean, Rcpp::Named("sd") = vec_sd);
}

double OrdinaryKriging::get_nllh(const Eigen::VectorXd& log_lengthscale, const double& nugget, double& mu, double& nu2, Eigen::LLT<Eigen::MatrixXd>& L, Eigen::VectorXd& a, Eigen::VectorXd& b) {
    Eigen::MatrixXd R = Ker_.compute(X_, log_lengthscale);
    R.diagonal().array() += nugget;
    L = R.llt();
    if (L.info() != 0) return 1e6;
    a = L.matrixL().solve(Eigen::VectorXd::Ones(n_));
    b = L.matrixL().solve(y_);
    mu = a.dot(b) / a.dot(a);
    nu2 = std::max(1.0 / (n_-1) * (b.dot(b) - mu*mu*a.dot(a)), 1e-15);
    double nllh = (n_-1) * std::log(nu2) + 2 * L.matrixLLT().diagonal().array().log().sum() + std::log(a.dot(a));
    if (interpolation_) {
        // penalty for non-interpolating
        Eigen::VectorXd ym = b - mu * a;
        Eigen::VectorXd pred = (mu + (L.matrixL()*ym - L.matrixU().solve(ym)*nugget).array()).matrix();
        double penalty = 2 * (y_-pred).squaredNorm() / (epsilon_ * y_tss_);
        nllh = nllh / n_ + penalty;
    }
    return nllh;
}

double OrdinaryKriging::get_nllh(const Eigen::VectorXd& log_lengthscale, const double& nugget) {
    double mu, nu2;
    Eigen::LLT<Eigen::MatrixXd> L(n_);
    Eigen::VectorXd a(n_), b(n_);
    double nllh = get_nllh(log_lengthscale, nugget, mu, nu2, L, a, b);
    return nllh;
}

void OrdinaryKriging::set_kriging_parameters() {
    nllh_ = get_nllh(Ker_.get_log_lengthscale(), nugget_, mu_, nu2_, L_, a_, b_);
}

void OrdinaryKriging::predict(const Eigen::VectorXd& xnew, double& mean, double& sd) {
    Eigen::VectorXd r = Ker_.compute(xnew, X_);
    Eigen::VectorXd d = L_.matrixL().solve(r);
    mean = mu_ + d.dot(b_-mu_*a_);
    sd = std::sqrt(nu2_) * std::sqrt(std::max(1 - d.dot(d) + std::pow(1-d.dot(a_),2)/a_.dot(a_), 0.0));
}

void LimitKriging::predict(const Eigen::VectorXd& xnew, double& mean, double& sd) {
    Eigen::VectorXd r = Ker_.compute(xnew, X_);
    Eigen::VectorXd d = L_.matrixL().solve(r);
    mean = d.dot(b_) / d.dot(a_);
    sd = std::sqrt(nu2_) * std::sqrt(std::max(1 - d.dot(d) + d.dot(d)*std::pow(1-d.dot(a_),2)/std::pow(d.dot(a_),2), 0.0));
}

RationalKriging::RationalKriging(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, Kernel& Ker, const bool& interpolation): Kriging(X,y,Ker,interpolation) {
    c_.resize(n_);
    s_.resize(n_);
}

void RationalKriging::add_data(const Eigen::MatrixXd& Xn, const Eigen::VectorXd& yn) {
    Kriging::add_data(Xn, yn);
    c_.resize(n_);
    s_.resize(n_);
}

double RationalKriging::get_nllh(const Eigen::VectorXd& log_lengthscale, const double& nugget, double& mu, double& nu2, Eigen::LLT<Eigen::MatrixXd>& L, Eigen::VectorXd& a, Eigen::VectorXd& b, Eigen::VectorXd& c, Eigen::VectorXd& s) {
    // compute covariance matrix 
    Eigen::MatrixXd R = Ker_.compute(X_, log_lengthscale);
    R.diagonal().array() += nugget;
    // get eigen vectors / values (rank from smallest to largest)
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigSolver(R);
    // compute gamma and c via binary search
    Eigen::VectorXd eigOne = eigSolver.eigenvectors().transpose() * Eigen::VectorXd::Ones(n_);
    double delta = eigSolver.eigenvalues().maxCoeff()/n_;
    double gamma, gamma_lb = 0.0, gamma_ub = 1.0, eps = 1e-4;
    c = eigSolver.eigenvectors() * (eigOne.array()/((1-gamma_lb)*eigSolver.eigenvalues().array()+gamma_lb)).matrix();
    if (c.minCoeff() > delta) {
        gamma = gamma_lb;
    } else {
        while (gamma_ub - gamma_lb > eps) {
            gamma = 0.5 * (gamma_lb + gamma_ub);
            c = eigSolver.eigenvectors() * (eigOne.array()/((1-gamma)*eigSolver.eigenvalues().array()+gamma)).matrix();
            if (c.minCoeff() < delta) {
                gamma_lb = gamma;
            } else {
                gamma_ub = gamma;
            }
        }
        gamma = gamma_ub;
        c = eigSolver.eigenvectors() * (eigOne.array()/((1-gamma)*eigSolver.eigenvalues().array()+gamma)).matrix();
    }
    c = c.cwiseMax(0).normalized();
    // compute negative likelihood
    s = R * c;
    L = R.llt();
    if (L.info() != 0) return 1e6;
    a = L.matrixL().solve(s);
    b = L.matrixL().solve(s.cwiseProduct(y_));
    mu = a.dot(b) / a.dot(a);
    nu2 = std::max(1.0 / (n_-1) * (b-mu*a).squaredNorm(), 1e-15);
    double nllh = (n_-1)*std::log(nu2) + 2*L.matrixLLT().diagonal().array().log().sum() - 2*s.array().log().sum() + std::log(c.dot(s));
    if (interpolation_) {
        // penalty for non-interpolating
        Eigen::VectorXd pred = (L.matrixL()*b - L.matrixU().solve(b)*nugget).cwiseQuotient(s);
        double penalty = 2 * (y_-pred).squaredNorm() / (epsilon_ * y_tss_);
        nllh = nllh / n_ + penalty;
    }
    return nllh;
}

double RationalKriging::get_nllh(const Eigen::VectorXd& log_lengthscale, const double& nugget) {
    double mu, nu2;
    Eigen::LLT<Eigen::MatrixXd> L(n_);
    Eigen::VectorXd a(n_), b(n_), c(n_), s(n_);
    double nllh = get_nllh(log_lengthscale, nugget, mu, nu2, L, a, b, c, s);
    return nllh;
}

void RationalKriging::set_kriging_parameters() {
    nllh_ = get_nllh(Ker_.get_log_lengthscale(), nugget_, mu_, nu2_, L_, a_, b_, c_, s_);
}

Eigen::VectorXd RationalKriging::get_c() {
    return c_;
} 

void RationalKriging::predict(const Eigen::VectorXd& xnew, double& mean, double& sd) {
    Eigen::VectorXd r = Ker_.compute(xnew, X_);
    Eigen::VectorXd d = L_.matrixL().solve(r);
    mean = d.dot(b_) / d.dot(a_);
    sd = std::sqrt(nu2_) * std::sqrt(std::max(1 - d.dot(d), 0.0)) / d.dot(a_);
}

GeneralizedRationalKriging::GeneralizedRationalKriging(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, Kernel& Ker, const bool& interpolation): Kriging(X,y,Ker,interpolation) {
    c_.resize(n_+1);
    s_.resize(n_);
}

void GeneralizedRationalKriging::add_data(const Eigen::MatrixXd& Xn, const Eigen::VectorXd& yn) {
    Kriging::add_data(Xn, yn);
    c_.resize(n_+1);
    s_.resize(n_);
}

double GeneralizedRationalKriging::get_nllh(const Eigen::VectorXd& log_lengthscale, const double& nugget, double& mu, double& nu2, Eigen::LLT<Eigen::MatrixXd>& L, Eigen::VectorXd& a, Eigen::VectorXd& b, Eigen::VectorXd& c, Eigen::VectorXd& s) {
    // compute covariance matrix 
    Eigen::MatrixXd R = Ker_.compute(X_, log_lengthscale);
    R.diagonal().array() += nugget;
    L = R.llt();
    if (L.info() != 0) return 1e6;
    // compute (c0, c')
    Eigen::MatrixXd M(n_+1, n_+1);
    M(0,0) = L.matrixL().solve(Eigen::VectorXd::Ones(n_)).squaredNorm();
    M.bottomLeftCorner(n_,1) = Eigen::VectorXd::Ones(n_);
    M.topRightCorner(1,n_) = Eigen::VectorXd::Ones(n_);
    M.bottomRightCorner(n_,n_) = R;
    Spectra::DenseSymMatProd<double> op(M);
    Spectra::SymEigsSolver<Spectra::DenseSymMatProd<double>> eigs(op, 1, (n_+1 < 20) ? n_+1 : 20);
    eigs.init();
    int nconv = eigs.compute(Spectra::SortRule::LargestAlge);
    c = eigs.eigenvectors().col(0).cwiseAbs();
    s = (c(0) + (R * c.tail(n_)).array()).matrix();
    a = L.matrixL().solve(s);
    b = L.matrixL().solve(s.cwiseProduct(y_));
    mu = a.dot(b) / a.dot(a);
    nu2 = std::max(1.0 / (n_-1) * (b-mu*a).squaredNorm(), 1e-15);
    double nllh = (n_-1)*std::log(nu2) + 2*L.matrixLLT().diagonal().array().log().sum() - 2*s.array().log().sum() + std::log(a.dot(a));
    if (interpolation_) {
        // penalty for non-interpolating
        Eigen::VectorXd ym = b - mu * a;
        Eigen::VectorXd pred = (mu + (L.matrixL()*ym - L.matrixU().solve(ym)*nugget).cwiseQuotient(s).array()).matrix();
        double penalty = 2 * (y_-pred).squaredNorm() / (epsilon_ * y_tss_);
        nllh = nllh / n_ + penalty;
    }
    return nllh;
}

double GeneralizedRationalKriging::get_nllh(const Eigen::VectorXd& log_lengthscale, const double& nugget) {
    double mu, nu2;
    Eigen::LLT<Eigen::MatrixXd> L(n_);
    Eigen::VectorXd a(n_), b(n_), c(n_+1), s(n_);
    double nllh = get_nllh(log_lengthscale, nugget, mu, nu2, L, a, b, c, s);
    return nllh;
}

void GeneralizedRationalKriging::set_kriging_parameters() {
    nllh_ = get_nllh(Ker_.get_log_lengthscale(), nugget_, mu_, nu2_, L_, a_, b_, c_, s_);
}

double GeneralizedRationalKriging::get_c0() {
    return c_(0);
}

Eigen::VectorXd GeneralizedRationalKriging::get_c() {
    return c_.tail(n_);
}

void GeneralizedRationalKriging::predict(const Eigen::VectorXd& xnew, double& mean, double& sd) {
    Eigen::VectorXd r = Ker_.compute(xnew, X_);
    Eigen::VectorXd d = L_.matrixL().solve(r);
    double sx = c_(0) + r.dot(c_.tail(n_));
    double tx = sx - d.dot(a_);
    mean = (mu_*tx + d.dot(b_)) / sx;
    sd = std::sqrt(nu2_) * std::sqrt(std::max(1 - d.dot(d) + std::pow(tx,2)/a_.dot(a_), 0.0)) / sx;
}

UniversalKriging::UniversalKriging(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, Kernel& Ker, const bool& interpolation, const std::size_t& pb, Rcpp::Function bfunc): Kriging(X,y,Ker,interpolation), pb_(pb), bfunc_(bfunc), LQ_(n_,pb_) {
    beta_.resize(pb_);
    B_.resize(n_, pb_);
    for (std::size_t i = 0; i < n_; i++) {
        std::vector<double> basis = Rcpp::as<std::vector<double>>(bfunc_(X_.row(i)));
        B_.row(i) = Eigen::Map<Eigen::VectorXd>(&basis[0], pb_);
    }
    // LQ_ = Eigen::ColPivHouseholderQR<Eigen::MatrixXd>(n_, pb_);
}

void UniversalKriging::add_data(const Eigen::MatrixXd& Xn, const Eigen::VectorXd& yn) {
    std::size_t nn = Xn.rows();
    Eigen::MatrixXd Bn(nn, pb_);
    for (std::size_t i = 0; i < nn; i++) {
        std::vector<double> basis = Rcpp::as<std::vector<double>>(bfunc_(Xn.row(i)));
        Bn.row(i) = Eigen::Map<Eigen::VectorXd>(&basis[0], pb_);
    }
    Eigen::Map<Eigen::MatrixXd> Bo(B_.data(), n_, pb_);
    B_.resize(n_+nn, pb_);
    B_.topRows(n_) = Bo;
    B_.bottomRows(nn) = Bn;      
    Kriging::add_data(Xn, yn);
    LQ_ = Eigen::ColPivHouseholderQR<Eigen::MatrixXd>(n_, pb_);
}

double UniversalKriging::get_nllh(const Eigen::VectorXd& log_lengthscale, const double& nugget, double& nu2, Eigen::LLT<Eigen::MatrixXd>& L, Eigen::ColPivHouseholderQR<Eigen::MatrixXd>& LQ, Eigen::VectorXd& a, Eigen::VectorXd& b, Eigen::VectorXd& beta) {
    Eigen::MatrixXd R = Ker_.compute(X_, log_lengthscale);
    R.diagonal().array() += nugget;
    L = R.llt();
    if (L.info() != 0) return 1e6;
    LQ = Eigen::ColPivHouseholderQR<Eigen::MatrixXd>(L.matrixL().solve(B_));
    b = L.matrixL().solve(y_);
    beta = LQ.solve(b);
    Eigen::VectorXd mu = B_ * beta;
    a = L.matrixL().solve(mu);
    nu2 = std::max(1.0 / (n_-pb_) * (b-a).squaredNorm(), 1e-15);
    double nllh = (n_-pb_) * std::log(nu2) + 2 * L.matrixLLT().diagonal().array().log().sum() + 2 * LQ.logAbsDeterminant();
    // penalty for non-interpolating
    if (interpolation_) {
        Eigen::VectorXd ym = b - a;
        Eigen::VectorXd pred = mu + (L.matrixL()*ym - L.matrixU().solve(ym)*nugget);
        double penalty = 2 * (y_-pred).squaredNorm() / (epsilon_ * std::min(y_tss_, (y_-mu).squaredNorm()-std::pow((y_-mu).sum(),2)/n_));
        nllh = nllh / n_ + penalty;
    }
    return nllh;
}

double UniversalKriging::get_nllh(const Eigen::VectorXd& log_lengthscale, const double& nugget) {
    double nu2;
    Eigen::LLT<Eigen::MatrixXd> L(n_);
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> LQ(n_,pb_);
    Eigen::VectorXd a(n_), b(n_), beta(pb_);
    double nllh = get_nllh(log_lengthscale, nugget, nu2, L, LQ, a, b, beta);
    return nllh;
}

void UniversalKriging::set_kriging_parameters() {
    nllh_ = get_nllh(Ker_.get_log_lengthscale(), nugget_, nu2_, L_, LQ_, a_, b_, beta_);
}

Eigen::VectorXd UniversalKriging::get_beta() {
    return beta_;
}

void UniversalKriging::predict(const Eigen::VectorXd& xnew, double& mean, double& sd) {
    std::vector<double> basis = Rcpp::as<std::vector<double>>(bfunc_(xnew));
    Eigen::Map<Eigen::VectorXd> bm(&basis[0], pb_);
    Eigen::VectorXd r = Ker_.compute(xnew, X_);
    Eigen::VectorXd d = L_.matrixL().solve(r);
    Eigen::VectorXd h = LQ_.colsPermutation().transpose()*bm - LQ_.matrixQR().triangularView<Eigen::Upper>().transpose()*(LQ_.householderQ().transpose()*d);
    mean = bm.dot(beta_) + d.dot(b_-a_);
    sd = std::sqrt(nu2_) * std::sqrt(std::max(1 - d.dot(d) + (LQ_.matrixQR().topLeftCorner(pb_,pb_).triangularView<Eigen::Upper>().transpose().solve(h)).squaredNorm(), 0.0));
}
