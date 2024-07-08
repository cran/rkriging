#include <limits>
#include <cmath>
#include <numeric>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include <boost/math/special_functions/bessel.hpp>
// [[Rcpp::depends(BH)]]
#include "kernel.h"

std::size_t Kernel::get_dimension() {
    return p_;
}

std::size_t Kernel::get_log_lengthscale_dimension() {
    return k_;
}

Eigen::VectorXd Kernel::get_lengthscale() {
    return get_lengthscale(get_log_lengthscale());
}

Eigen::MatrixXd Kernel::compute(const Eigen::MatrixXd& X) {
    return compute(X, get_log_lengthscale());
}

Eigen::MatrixXd Kernel::compute(const Eigen::MatrixXd& X, const Eigen::VectorXd& log_lengthscale) {
    Eigen::VectorXd lengthscale = get_lengthscale(log_lengthscale);
    std::size_t n = X.rows(), p = X.cols();
    Eigen::MatrixXd R(n,n);
    Eigen::VectorXd sqdist(p);
    for (std::size_t i = 0; i < n; i++) {
        for (std::size_t j = i; j < n; j++) {
            if (i == j) {
                R(i,j) = 1.0;
            } else {
                sqdist = (X.row(i)-X.row(j)).transpose().cwiseQuotient(lengthscale).cwiseAbs2();
                R(i,j) = evaluate(sqdist);
                R(j,i) = R(i,j);
            }
        }
    }
    return R;
}

Eigen::VectorXd Kernel::compute(const Eigen::VectorXd& x, const Eigen::MatrixXd& X) {
    Eigen::VectorXd lengthscale = get_lengthscale();
    std::size_t n = X.rows(), p = X.cols();
    Eigen::VectorXd r(n);
    Eigen::VectorXd sqdist(p);
    for (std::size_t i = 0; i < n; i++) {
        sqdist = (x-X.row(i).transpose()).cwiseQuotient(lengthscale).cwiseAbs2();
        r(i) = evaluate(sqdist);
    }
    return r;
}

BaseKernel::BaseKernel(const Eigen::VectorXd& lengthscale) {
    p_ = lengthscale.size();
    k_ = lengthscale.size();
    log_lengthscale_.resize(k_);
    set_lengthscale(lengthscale);
}

BaseKernel::BaseKernel(const double& lengthscale, const Eigen::VectorXd& weight): weight_(weight) {
    p_ = weight_.size();
    k_ = 1;
    log_lengthscale_.resize(k_);
    log_lengthscale_[0] = std::log(lengthscale);
}

void BaseKernel::set_log_lengthscale(const Eigen::VectorXd& log_lengthscale) {
    log_lengthscale_ = log_lengthscale;
}

Eigen::VectorXd BaseKernel::get_log_lengthscale() {
    return log_lengthscale_;
}

void BaseKernel::set_lengthscale(const Eigen::VectorXd& lengthscale) {
    if (k_ < p_) {
        Eigen::VectorXd w = lengthscale.array().square().inverse().matrix();
        double wsum = w.sum();
        log_lengthscale_[0] = -0.5 * std::log(wsum);
        weight_ = w / wsum;
    } else {
        log_lengthscale_ = lengthscale.array().log().matrix();
    }
}

Eigen::VectorXd BaseKernel::get_lengthscale(const Eigen::VectorXd& log_lengthscale) {
    if (k_ < p_) {
        return std::exp(log_lengthscale[0]) * weight_.cwiseSqrt();
    } else {
        return log_lengthscale.array().exp().matrix();
    }
}

void BaseKernel::relax_lengthscale_constraint() {
    Eigen::VectorXd lengthscale = get_lengthscale(get_log_lengthscale());
    k_ = p_;
    log_lengthscale_.resize(k_);
    set_lengthscale(lengthscale);
}

double BaseKernel::evaluate(const Eigen::VectorXd& sqdist) {
    return evaluate(sqdist.sum());
}

double GaussianKernel::evaluate(const double& sqdist) {
    return std::exp(-0.5*sqdist);
}

double RQKernel::evaluate(const double& sqdist) {
    return std::pow(1.0+sqdist/(2.0*alpha_),-alpha_);
}

double Matern12Kernel::evaluate(const double& sqdist) {
    return std::exp(-std::sqrt(sqdist));
}

double Matern32Kernel::evaluate(const double& sqdist) {
    return (1+std::sqrt(3*sqdist)) * std::exp(-std::sqrt(3*sqdist));
}

double Matern52Kernel::evaluate(const double& sqdist) {
    return (1+std::sqrt(5*sqdist)+5*sqdist/3) * std::exp(-std::sqrt(5*sqdist));   
}

double MaternKernel::evaluate(const double& sqdist) {
    double scaled_dist = std::sqrt(2*nu_) * std::sqrt(sqdist);
    scaled_dist = std::max(1e-10, scaled_dist);
    return std::pow(2.0,1-nu_) / std::tgamma(nu_) * std::pow(scaled_dist,nu_) * boost::math::cyl_bessel_k(nu_,scaled_dist);
}

double UDFKernel::evaluate(const double& sqdist) {
    Rcpp::NumericVector res = kfunc_(sqdist);
    return res(0);
}

double MultiplicativeRQKernel::evaluate(const Eigen::VectorXd& sqdist) {
    return sqdist.unaryExpr([&](double x){return evaluate(x);}).prod();
}

double MultiplicativeMaternKernel::evaluate(const Eigen::VectorXd& sqdist) {
    return sqdist.unaryExpr([&](double x){return evaluate(x);}).prod();
}

double MultiplicativeUDFKernel::evaluate(const Eigen::VectorXd& sqdist) {
    return sqdist.unaryExpr([&](double x){return evaluate(x);}).prod();
}
