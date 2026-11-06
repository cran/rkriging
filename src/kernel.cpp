#include <limits>
#include <cmath>
#include <numeric>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include <boost/math/special_functions/bessel.hpp>
// [[Rcpp::depends(BH)]]
#include "kernel.h"

std::size_t Kernel::get_dimension() const {
    return _p;
}

std::size_t Kernel::get_log_lengthscale_dimension() const {
    return _k;
}

Eigen::VectorXd Kernel::get_lengthscale() const {
    return _get_lengthscale(get_log_lengthscale());
}

Eigen::MatrixXd Kernel::compute(const Eigen::MatrixXd& X) const {
    return compute(X, get_log_lengthscale());
}

Eigen::MatrixXd Kernel::compute(const Eigen::MatrixXd& X, const Eigen::VectorXd& log_lengthscale) const {
    Eigen::VectorXd lengthscale = _get_lengthscale(log_lengthscale);
    std::size_t n = X.rows(), p = X.cols();
    Eigen::MatrixXd R(n, n);
    Eigen::VectorXd sqdist(p);
    for (std::size_t i = 0; i < n; i++) {
        for (std::size_t j = i; j < n; j++) {
            if (i == j) {
                R(i, j) = 1.0;
            } else {
                sqdist = (X.row(i) - X.row(j)).transpose().cwiseQuotient(lengthscale).cwiseAbs2();
                R(i, j) = _evaluate(sqdist);
                R(j, i) = R(i, j);
            }
        }
    }
    return R;
}

Eigen::VectorXd Kernel::compute(const Eigen::VectorXd& x, const Eigen::MatrixXd& X) const {
    Eigen::VectorXd lengthscale = _get_lengthscale(get_log_lengthscale());
    std::size_t n = X.rows(), p = X.cols();
    Eigen::VectorXd r(n);
    Eigen::VectorXd sqdist(p);
    for (std::size_t i = 0; i < n; i++) {
        sqdist = (x - X.row(i).transpose()).cwiseQuotient(lengthscale).cwiseAbs2();
        r(i) = _evaluate(sqdist);
    }
    return r;
}

BaseKernel::BaseKernel(const Eigen::VectorXd& lengthscale) {
    _p = lengthscale.size();
    _k = lengthscale.size();
    _log_lengthscale.resize(_k);
    set_lengthscale(lengthscale);
}

BaseKernel::BaseKernel(const double& lengthscale, const Eigen::VectorXd& weight): _weight(weight) {
    _p = _weight.size();
    _k = 1;
    _log_lengthscale.resize(_k);
    _log_lengthscale[0] = std::log(lengthscale);
}

void BaseKernel::set_log_lengthscale(const Eigen::VectorXd& log_lengthscale) {
    _log_lengthscale = log_lengthscale;
}

Eigen::VectorXd BaseKernel::get_log_lengthscale() const {
    return _log_lengthscale;
}

void BaseKernel::set_lengthscale(const Eigen::VectorXd& lengthscale) {
    if (lengthscale.size() != _k)
        Rcpp::stop("Lengthscale dimension does not match kernel dimension.");
    if ((lengthscale.array() <= 0).any())
        Rcpp::stop("Lengthscale must be positive.");
    if (_k < _p) {
        Eigen::VectorXd w = lengthscale.array().square().inverse().matrix();
        double wsum = w.sum();
        _log_lengthscale[0] = -0.5 * std::log(wsum);
        _weight = w / wsum;
    } else {
        _log_lengthscale = lengthscale.array().log().matrix();
    }
}

Eigen::VectorXd BaseKernel::_get_lengthscale(const Eigen::VectorXd& log_lengthscale) const {
    if (_k < _p) {
        return std::exp(log_lengthscale[0]) * _weight.cwiseSqrt();
    } else {
        return log_lengthscale.array().exp().matrix();
    }
}

void BaseKernel::relax_lengthscale_constraint() {
    Eigen::VectorXd lengthscale = _get_lengthscale(get_log_lengthscale());
    _k = _p;
    _log_lengthscale.resize(_k);
    set_lengthscale(lengthscale);
}

double BaseKernel::_evaluate(const Eigen::VectorXd& sqdist) const {
    return _evaluate(sqdist.sum());
}

double GaussianKernel::_evaluate(const double& sqdist) const {
    return std::exp(-0.5*sqdist);
}

double RQKernel::_evaluate(const double& sqdist) const {
    return std::pow(1.0 + sqdist / (2.0 * _alpha), -_alpha);
}

double IMKernel::_evaluate(const double& sqdist) const {
    return std::pow(1.0 + sqdist / 2.0, -_alpha);
}

double Matern12Kernel::_evaluate(const double& sqdist) const {
    return std::exp(-std::sqrt(sqdist));
}

double Matern32Kernel::_evaluate(const double& sqdist) const {
    return (1 + std::sqrt(3 * sqdist)) * std::exp(-std::sqrt(3 * sqdist));
}

double Matern52Kernel::_evaluate(const double& sqdist) const {
    return (1 + std::sqrt(5 * sqdist) + 5 * sqdist / 3) * std::exp(-std::sqrt(5 * sqdist));
}

double MaternKernel::_evaluate(const double& sqdist) const {
    double scaled_dist = std::sqrt(2 * _nu) * std::sqrt(sqdist);
    scaled_dist = std::max(1e-10, scaled_dist);
    return std::pow(2.0, 1 - _nu) / std::tgamma(_nu) 
         * std::pow(scaled_dist, _nu) 
         * boost::math::cyl_bessel_k(_nu, scaled_dist);
}

double UDFKernel::_evaluate(const double& sqdist) const {
    Rcpp::NumericVector res = _kfunc(sqdist);
    return res(0);
}
