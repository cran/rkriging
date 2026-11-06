#ifndef KERNEL_H
#define KERNEL_H

#include <limits>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

class Kernel {
    protected:
        std::size_t _p;
        std::size_t _k;
        virtual double _evaluate(const Eigen::VectorXd& sqdist) const = 0;
        virtual Eigen::VectorXd _get_lengthscale(const Eigen::VectorXd& log_lengthscale) const = 0;
    public:
        Kernel() = default;
        std::size_t get_dimension() const;
        std::size_t get_log_lengthscale_dimension() const;
        Eigen::VectorXd get_lengthscale() const;
        virtual void set_log_lengthscale(const Eigen::VectorXd& log_lengthscale) = 0;
        virtual Eigen::VectorXd get_log_lengthscale() const = 0;
        virtual void set_lengthscale(const Eigen::VectorXd& lengthscale) = 0;
        virtual void relax_lengthscale_constraint() = 0;
        Eigen::MatrixXd compute(const Eigen::MatrixXd& X) const;
        Eigen::MatrixXd compute(const Eigen::MatrixXd& X, const Eigen::VectorXd& log_lengthscale) const;
        Eigen::VectorXd compute(const Eigen::VectorXd& x, const Eigen::MatrixXd& X) const;
        virtual ~Kernel() = default;
};

class BaseKernel : public Kernel {
    protected:
        Eigen::VectorXd _log_lengthscale;
        Eigen::VectorXd _weight;
        double _evaluate(const Eigen::VectorXd& sqdist) const override;
        Eigen::VectorXd _get_lengthscale(const Eigen::VectorXd& log_lengthscale) const override;
        virtual double _evaluate(const double& sqdist) const = 0;
    public:
        BaseKernel(const Eigen::VectorXd& lengthscale);
        BaseKernel(const double& lengthscale, const Eigen::VectorXd& weight);
        void set_log_lengthscale(const Eigen::VectorXd& log_lengthscale) override;
        Eigen::VectorXd get_log_lengthscale() const override;
        void set_lengthscale(const Eigen::VectorXd& lengthscale) override;
        void relax_lengthscale_constraint() override;
};

class GaussianKernel : public BaseKernel {
    protected:
        double _evaluate(const double& sqdist) const override;
    public:
        using BaseKernel::BaseKernel;
};

class RQKernel : public BaseKernel {
    protected:
        double _alpha;
        double _evaluate(const double& sqdist) const override;
    public:
        RQKernel(const Eigen::VectorXd& lengthscale, const double& alpha)
            : BaseKernel(lengthscale), _alpha(alpha) {};
        RQKernel(const double& lengthscale, const Eigen::VectorXd& weight, const double& alpha)
            : BaseKernel(lengthscale, weight), _alpha(alpha) {};
};

class IMKernel : public RQKernel {
    protected:
        double _evaluate(const double& sqdist) const override;
    public:
        using RQKernel::RQKernel;
};

class Matern12Kernel : public BaseKernel {
    protected:
        double _evaluate(const double& sqdist) const override;
    public:
        using BaseKernel::BaseKernel;
};

class Matern32Kernel : public BaseKernel {
    protected:
        double _evaluate(const double& sqdist) const override;
    public:
        using BaseKernel::BaseKernel;
};

class Matern52Kernel : public BaseKernel {
    protected:
        double _evaluate(const double& sqdist) const override;
    public:
        using BaseKernel::BaseKernel;
};

class MaternKernel : public BaseKernel {
    protected:
        double _nu;
        double _evaluate(const double& sqdist) const override;
    public:
        MaternKernel(const Eigen::VectorXd& lengthscale, const double& nu)
            : BaseKernel(lengthscale), _nu(nu) {};
        MaternKernel(const double& lengthscale, const Eigen::VectorXd& weight, const double& nu)
            : BaseKernel(lengthscale, weight), _nu(nu) {};
};

class UDFKernel : public BaseKernel {
    protected:
        Rcpp::Function _kfunc;
        double _evaluate(const double& sqdist) const override;
    public:
        UDFKernel(const Eigen::VectorXd& lengthscale, Rcpp::Function kfunc)
            : BaseKernel(lengthscale), _kfunc(kfunc) {};
        UDFKernel(const double& lengthscale, const Eigen::VectorXd& weight, Rcpp::Function kfunc)
            : BaseKernel(lengthscale, weight), _kfunc(kfunc) {};
};

template<typename BaseKernelType>
class MultiplicativeKernel : public BaseKernelType {
    protected:
        using BaseKernelType::_evaluate;
        double _evaluate(const Eigen::VectorXd& sqdist) const override {
            return sqdist.unaryExpr([this](double x) {
                return BaseKernelType::_evaluate(x);
            }).prod();
        }
    public:
        using BaseKernelType::BaseKernelType;
};

using MultiplicativeRQKernel = MultiplicativeKernel<RQKernel>;
using MultiplicativeIMKernel = MultiplicativeKernel<IMKernel>;
using MultiplicativeMaternKernel = MultiplicativeKernel<MaternKernel>;
using MultiplicativeUDFKernel = MultiplicativeKernel<UDFKernel>;

#endif
