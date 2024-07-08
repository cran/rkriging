#ifndef KERNEL_H
#define KERNEL_H

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

class Kernel {
    protected:
        std::size_t p_;
        std::size_t k_;
    public:
        Kernel() {};
        std::size_t get_dimension();
        std::size_t get_log_lengthscale_dimension();
        Eigen::VectorXd get_lengthscale();
        virtual void set_log_lengthscale(const Eigen::VectorXd& log_lengthscale) = 0;
        virtual Eigen::VectorXd get_log_lengthscale() = 0;
        virtual void set_lengthscale(const Eigen::VectorXd& lengthscale) = 0;
        virtual Eigen::VectorXd get_lengthscale(const Eigen::VectorXd& log_lengthscale) = 0;
        virtual void relax_lengthscale_constraint() = 0;
        virtual double evaluate(const Eigen::VectorXd& sqdist) = 0;
        Eigen::MatrixXd compute(const Eigen::MatrixXd& X);
        Eigen::MatrixXd compute(const Eigen::MatrixXd& X, const Eigen::VectorXd& log_lengthscale);
        Eigen::VectorXd compute(const Eigen::VectorXd& x, const Eigen::MatrixXd& X);
        ~Kernel() {};
};

class BaseKernel : public Kernel {
    protected:
        Eigen::VectorXd log_lengthscale_;
        Eigen::VectorXd weight_;
    public:
        BaseKernel(const Eigen::VectorXd& lengthscale);
        BaseKernel(const double& lengthscale, const Eigen::VectorXd& weight);
        void set_log_lengthscale(const Eigen::VectorXd& log_lengthscale) override;
        Eigen::VectorXd get_log_lengthscale() override;
        void set_lengthscale(const Eigen::VectorXd& lengthscale) override;
        Eigen::VectorXd get_lengthscale(const Eigen::VectorXd& log_lengthscale) override;
        void relax_lengthscale_constraint() override;
        double evaluate(const Eigen::VectorXd& sqdist) override;
        virtual double evaluate(const double& sqdist) = 0;
};

class GaussianKernel : public BaseKernel {
    public:
        using BaseKernel::BaseKernel;
        double evaluate(const double& sqdist) override;
 };

class RQKernel : public BaseKernel {
    protected:
        double alpha_;
    public:
        RQKernel(const Eigen::VectorXd& lengthscale, const double& alpha): BaseKernel(lengthscale), alpha_(alpha) {};
        RQKernel(const double& lengthscale, const Eigen::VectorXd& weight, const double& alpha): BaseKernel(lengthscale,weight), alpha_(alpha) {};
        double evaluate(const double& sqdist) override;
};

class Matern12Kernel : public BaseKernel {
    public:
        using BaseKernel::BaseKernel;
        double evaluate(const double& sqdist) override;
};

class Matern32Kernel : public BaseKernel {
    public:
        using BaseKernel::BaseKernel;
        double evaluate(const double& sqdist) override;
};

class Matern52Kernel : public BaseKernel {
    public:
        using BaseKernel::BaseKernel;
        double evaluate(const double& sqdist) override;
};

class MaternKernel : public BaseKernel {
    protected:
        double nu_;
    public:
        MaternKernel(const Eigen::VectorXd& lengthscale, const double& nu): BaseKernel(lengthscale), nu_(nu) {};
        MaternKernel(const double& lengthscale, const Eigen::VectorXd& weight, const double& nu): BaseKernel(lengthscale,weight), nu_(nu) {};
        double evaluate(const double& sqdist) override;
};

class UDFKernel : public BaseKernel {
    protected:
        Rcpp::Function kfunc_;
    public:
        UDFKernel(const Eigen::VectorXd& lengthscale, Rcpp::Function kfunc): BaseKernel(lengthscale), kfunc_(kfunc) {};
        UDFKernel(const double& lengthscale, const Eigen::VectorXd& weight, Rcpp::Function kfunc): BaseKernel(lengthscale,weight), kfunc_(kfunc) {};
        double evaluate(const double& sqdist) override;
};

class MultiplicativeRQKernel : public RQKernel {
    public:
        using RQKernel::RQKernel;
        using RQKernel::evaluate;
        double evaluate(const Eigen::VectorXd& sqdist) override;
};

class MultiplicativeMaternKernel : public MaternKernel {
    public:
        using MaternKernel::MaternKernel;
        using MaternKernel::evaluate;
        double evaluate(const Eigen::VectorXd& sqdist) override;
};

class MultiplicativeUDFKernel : public UDFKernel {
    public:
        using UDFKernel::UDFKernel;
        using UDFKernel::evaluate;
        double evaluate(const Eigen::VectorXd& sqdist) override;
};

#endif
