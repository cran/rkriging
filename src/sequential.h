#ifndef SEQUENTIAL_H
#define SEQUENTIAL_H

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include "utils.h"
#include "kriging.h"

double _nlopt_criterion(unsigned n, const double* x, double* grad, void* data);

class Sequential {
    protected:
        Kriging& _Krig;
        Eigen::VectorXd _design_lower_bound;
        Eigen::VectorXd _design_upper_bound;
        std::string _nlopt_algorithm = "NLOPT_GN_MLSL_LDS";
        std::string _nlopt_local_algorithm = "NLOPT_LN_SBPLX";
        std::size_t _nlopt_maxeval = 1000;
        std::size_t _nlopt_population = 10;
        virtual void _add_data(const Eigen::MatrixXd& X, const Eigen::VectorXd& y);
        virtual double _get_criterion(const double& mean, const double& sd) = 0;
    public:
        Sequential(Kriging& Krig,
                   const Eigen::VectorXd& design_lower_bound,
                   const Eigen::VectorXd& design_upper_bound);
        void add_data(const Eigen::MatrixXd& X, const Eigen::VectorXd& y);
        void fit_kriging(bool optimize);
        double get_criterion(const Eigen::VectorXd& design);
        Rcpp::List get_data() const;
        std::size_t get_design_dimension() const;
        Eigen::VectorXd get_design_lower_bound() const;
        Eigen::VectorXd get_design_upper_bound() const;
        Eigen::VectorXd get_next_design(Eigen::VectorXd& x);
        Eigen::VectorXd get_next_design_from_candidates(const Eigen::MatrixXd& Xc);
        Eigen::VectorXd predict(const Eigen::VectorXd& xnew) const;
        Rcpp::List predict(const Eigen::MatrixXd& Xnew) const;
        void set_design_bound(const Eigen::VectorXd& design_lower_bound,
                              const Eigen::VectorXd& design_upper_bound);
        void set_nlopt_algorithm(const std::string& nlopt_algorithm);
        void set_nlopt_local_algorithm(const std::string& nlopt_local_algorithm);
        void set_nlopt_maxeval(const std::size_t& nlopt_maxeval);
        void set_nlopt_population(const std::size_t& nlopt_population);
        virtual ~Sequential() = default;
};

class EmulationSequential : public Sequential {
    public:
        using Sequential::Sequential;
};

class OptimizationSequential : public Sequential {
    protected:
        double _ymin;
        std::size_t _ymin_index;
        void _add_data(const Eigen::MatrixXd& X, const Eigen::VectorXd& y) override;
    public:
        OptimizationSequential(Kriging& Krig,
                               const Eigen::VectorXd& design_lower_bound,
                               const Eigen::VectorXd& design_upper_bound);
        double get_best_objective() const;
        Eigen::VectorXd get_best_design() const;
};

class ALMSequential : public EmulationSequential {
    protected:
        double _get_criterion(const double& mean, const double& sd) override;
    public:
        using EmulationSequential::EmulationSequential;
};

class UCBSequential : public OptimizationSequential {
    protected:
        double _lambda;
        double _get_criterion(const double& mean, const double& sd) override;
    public:
        UCBSequential(Kriging& Krig,
                      const Eigen::VectorXd& design_lower_bound,
                      const Eigen::VectorXd& design_upper_bound,
                      const double& lambda)
            : OptimizationSequential(Krig, design_lower_bound, design_upper_bound),
              _lambda(lambda) {};
};

class PISequential : public OptimizationSequential {
    protected:
        double _get_criterion(const double& mean, const double& sd) override;
    public:
        using OptimizationSequential::OptimizationSequential;
};

class EISequential : public OptimizationSequential {
    protected:
        double _get_criterion(const double& mean, const double& sd) override;
    public:
        using OptimizationSequential::OptimizationSequential;
};

class UDFSequential : public Sequential {
    protected:
        Rcpp::Function _cfunc;
        Rcpp::List _cparams;
        double _get_criterion(const double& mean, const double& sd) override;
    public:
        UDFSequential(Kriging& Krig,
                      const Eigen::VectorXd& design_lower_bound,
                      const Eigen::VectorXd& design_upper_bound,
                      Rcpp::Function cfunc, Rcpp::List cparams)
            : Sequential(Krig, design_lower_bound, design_upper_bound),
              _cfunc(cfunc), _cparams(cparams) {};
};

#endif
