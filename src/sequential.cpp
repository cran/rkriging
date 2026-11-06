#include <limits>
#include <cmath>
#include <vector>
#include <numeric>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include <nloptrAPI.h>
// [[Rcpp::depends(nloptr)]]
#include "utils.h"
#include "kriging.h"
#include "sequential.h"

double _nlopt_criterion(unsigned n, const double* x, double* grad, void* data) {
    Sequential* sequential = static_cast<Sequential*>(data);
    Eigen::Map<Eigen::VectorXd> design(const_cast<double*>(x), sequential->get_design_dimension());
    return sequential->get_criterion(design);
}

Sequential::Sequential(Kriging& Krig,
                       const Eigen::VectorXd& design_lower_bound,
                       const Eigen::VectorXd& design_upper_bound)
    : _Krig(Krig), _design_lower_bound(design_lower_bound),
      _design_upper_bound(design_upper_bound) {
    if (!Krig.is_interpolation())
        Rcpp::stop("Sequential design only supports interpolation kriging.");
}

std::size_t Sequential::get_design_dimension() const {
    return _Krig.get_dimension();
}

Eigen::VectorXd Sequential::get_design_lower_bound() const {
    return _design_lower_bound;
}

Eigen::VectorXd Sequential::get_design_upper_bound() const {
    return _design_upper_bound;
}

void Sequential::set_design_bound(const Eigen::VectorXd& design_lower_bound,
                                  const Eigen::VectorXd& design_upper_bound) {
    _design_lower_bound = design_lower_bound;
    _design_upper_bound = design_upper_bound;
}

void Sequential::set_nlopt_algorithm(const std::string& nlopt_algorithm) {
    _nlopt_algorithm = nlopt_algorithm;
}

void Sequential::set_nlopt_local_algorithm(const std::string& nlopt_local_algorithm) {
    _nlopt_local_algorithm = nlopt_local_algorithm;
}

void Sequential::set_nlopt_maxeval(const std::size_t& nlopt_maxeval) {
    _nlopt_maxeval = nlopt_maxeval;
}

void Sequential::set_nlopt_population(const std::size_t& nlopt_population) {
    _nlopt_population = nlopt_population;
}

double Sequential::get_criterion(const Eigen::VectorXd& design) {
    Eigen::VectorXd pred = _Krig.predict(design);
    return _get_criterion(pred(0), pred(1));
}

Eigen::VectorXd Sequential::get_next_design(Eigen::VectorXd& x) {
    if (x.size() != get_design_dimension())
        Rcpp::stop("Starting design dimension does not match kriging dimension.");
    if ((x.array() < _design_lower_bound.array()).any() || (x.array() > _design_upper_bound.array()).any())
        Rcpp::stop("Starting design is outside the design bounds.");
    std::size_t opt_dim = x.size();
    Eigen::Map<Eigen::VectorXd> lb(_design_lower_bound.data(), opt_dim);
    Eigen::Map<Eigen::VectorXd> ub(_design_upper_bound.data(), opt_dim);

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
    nlopt_set_max_objective(optimizer, _nlopt_criterion, this);
    nlopt_set_maxeval(optimizer, _nlopt_maxeval);

    double criterion_max;
    nlopt_optimize(optimizer, x.data(), &criterion_max);
    nlopt_destroy(optimizer);

    return x;
}

Eigen::VectorXd Sequential::get_next_design_from_candidates(const Eigen::MatrixXd& Xc) {
    if (Xc.cols() != get_design_dimension())
        Rcpp::stop("Candidate designs dimension does not match kriging dimension.");
    for (std::size_t j = 0; j < Xc.cols(); j++) {
        if ((Xc.col(j).array() < _design_lower_bound(j)).any() ||
            (Xc.col(j).array() > _design_upper_bound(j)).any())
            Rcpp::stop("Candidate designs are outside the design bounds.");
    }
    double criterion, max_criterion = 0;
    std::size_t max_index = 0;
    for (std::size_t i = 0; i < Xc.rows(); i++) {
        criterion = get_criterion(Xc.row(i));
        if (criterion > max_criterion) {
            max_criterion = criterion;
            max_index = i;
        }
    }
    return Xc.row(max_index);
}

void Sequential::fit_kriging(bool optimize) {
    _Krig.fit(optimize);
}

void Sequential::_add_data(const Eigen::MatrixXd& X, const Eigen::VectorXd& y) {
    _Krig.add_data(X, y);
}

void Sequential::add_data(const Eigen::MatrixXd& X, const Eigen::VectorXd& y) {
    _add_data(X, y);
}

Rcpp::List Sequential::get_data() const {
    return _Krig.get_data();
}

Eigen::VectorXd Sequential::predict(const Eigen::VectorXd& xnew) const {
    return _Krig.predict(xnew);
}

Rcpp::List Sequential::predict(const Eigen::MatrixXd& Xnew) const {
    return _Krig.predict(Xnew);
}

OptimizationSequential::OptimizationSequential(Kriging& Krig,
                                               const Eigen::VectorXd& design_lower_bound,
                                               const Eigen::VectorXd& design_upper_bound)
    : Sequential(Krig, design_lower_bound, design_upper_bound) {
    Eigen::Index index;
    _ymin = Krig.get_y().minCoeff(&index);
    _ymin_index = static_cast<std::size_t>(index);
}

double OptimizationSequential::get_best_objective() const {
    return _ymin;
}

Eigen::VectorXd OptimizationSequential::get_best_design() const {
    return _Krig.get_X().row(_ymin_index);
}

void OptimizationSequential::_add_data(const Eigen::MatrixXd& X, const Eigen::VectorXd& y) {
    Eigen::Index index;
    double ymin = y.minCoeff(&index);
    if (ymin < _ymin) {
        _ymin = ymin;
        _ymin_index = _Krig.get_datasize() + static_cast<std::size_t>(index);
    }
    Sequential::_add_data(X, y);
}

double ALMSequential::_get_criterion(const double& mean, const double& sd) {
    return sd;
}

double UCBSequential::_get_criterion(const double& mean, const double& sd) {
    return -(mean - _lambda * sd);
}

double PISequential::_get_criterion(const double& mean, const double& sd) {
    double u = (_ymin - mean) / sd;
    return R::pnorm(u, 0, 1, true, false);
}

double EISequential::_get_criterion(const double& mean, const double& sd) {
    double u = (_ymin - mean) / sd;
    return sd * (u * R::pnorm(u, 0, 1, true, false) + R::dnorm(u, 0, 1, false));
}

double UDFSequential::_get_criterion(const double& mean, const double& sd) {
    Rcpp::NumericVector res = _cfunc(mean, sd, _Krig.get_X(), _Krig.get_y(), _cparams);
    return res(0);
}
