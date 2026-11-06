#ifndef KRINGING_H
#define KRINGING_H

#include <limits>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include "utils.h"
#include "kernel.h"

// Global constants for kriging
const double STABILITY_THRESHOLD = 1e-12;

double _nlopt_loss(unsigned n, const double* x, double* grad, void* data);

class Kriging {
    protected:
        std::size_t _n;
        std::size_t _p;
        Eigen::MatrixXd _X;
        Eigen::VectorXd _y;
        Kernel& _Ker;
        bool _interpolation;
        double _log_nugget = std::log(STABILITY_THRESHOLD);
        double _epsilon = 1e-3;
        double _mu;
        double _nu2;
        Eigen::VectorXd _a;
        Eigen::VectorXd _b;
        Eigen::LLT<Eigen::MatrixXd> _L;
        Eigen::MatrixXd _Kinv_pred;
        double _loss = std::numeric_limits<double>::infinity();
        double _y_tss;
        std::string _nlopt_algorithm = "NLOPT_LN_SBPLX";
        std::string _nlopt_local_algorithm = "";
        std::size_t _nlopt_maxeval = 100;
        std::size_t _nlopt_population = 10;
        std::string _loss_type = "NLLH";
        Eigen::VectorXd _log_lengthscale_lower_bound;
        Eigen::VectorXd _log_lengthscale_upper_bound;
        double _log_nugget_lower_bound = std::log(STABILITY_THRESHOLD);
        double _log_nugget_upper_bound = 0;
        bool _is_fitted = false;
        // Private helper methods
        void _invalidate_fit();
        void _raise_error_if_not_fitted() const;
        void _set_default_lengthscale_bound();
        void _set_y_tss();
        void _validate_lengthscale_bound();
        void _validate_nugget_bound();
        void _validate_loss_type(const std::string& loss_type);
        void _update_prediction_inverse_neumann(const std::size_t& max_iter = 5);
        // Virtual methods for subclass-specific logic
        virtual void _add_data(const Eigen::MatrixXd& Xn, const Eigen::VectorXd& yn);
        virtual double _get_nllh() const = 0;
        virtual double _get_loocv() const = 0;
        virtual void _predict(const Eigen::VectorXd& xnew, double& mean, double& sd) const = 0;
        virtual Rcpp::List _predict_neumann(const Eigen::MatrixXd& Xnew) const = 0;
        virtual bool _set_hyperparameters(const Eigen::VectorXd& log_lengthscale, const double& log_nugget) = 0;
        // Private internals called by fit()
        void _optimize_hyperparameters();
        void _set_hyperparameters();
    public:
        Kriging(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, Kernel& Ker,
                const bool& interpolation, const std::string& loss_type);
        void add_data(const Eigen::MatrixXd& Xn, const Eigen::VectorXd& yn);
        double get_mu() const;
        double get_nu2() const;
        double get_sigma2() const;
        double get_nugget() const;
        double get_log_nugget() const;
        Rcpp::List get_data() const;
        Eigen::MatrixXd get_X() const;
        Eigen::VectorXd get_y() const;
        std::size_t get_datasize() const;
        std::size_t get_dimension() const;
        Eigen::VectorXd get_lengthscale() const;
        std::size_t get_lengthscale_dimension() const;
        double get_loss() const;
        double get_loss(const Eigen::VectorXd& log_lengthscale, const double& log_nugget);
        std::string get_loss_type() const;
        double get_nllh() const;
        double get_loocv() const;
        bool is_fitted() const;
        bool is_interpolation() const;
        void relax_lengthscale_constraint();
        void set_epsilon(const double& epsilon);
        void set_lengthscale(const Eigen::VectorXd& lengthscale);
        void set_lengthscale_bound(const Eigen::VectorXd& lengthscale_lower_bound,
                                   const Eigen::VectorXd& lengthscale_upper_bound);
        void set_nugget(const double& nugget);
        void set_nugget_bound(const double& nugget_lower_bound,
                              const double& nugget_upper_bound);
        void set_nlopt_algorithm(const std::string& nlopt_algorithm);
        void set_nlopt_local_algorithm(const std::string& nlopt_local_algorithm);
        void set_nlopt_maxeval(const std::size_t& nlopt_maxeval);
        void set_nlopt_population(const std::size_t& nlopt_population);
        void fit(bool optimize);
        Eigen::VectorXd predict(const Eigen::VectorXd& xnew) const;
        Rcpp::List predict(const Eigen::MatrixXd& Xnew) const;
        Rcpp::List predict_neumann(const Eigen::MatrixXd& Xnew) const;
        virtual ~Kriging() = default;
};

class OrdinaryKriging : public Kriging {
    protected:
        double _get_nllh() const override;
        double _get_loocv() const override;
        void _predict(const Eigen::VectorXd& xnew, double& mean, double& sd) const override;
        Rcpp::List _predict_neumann(const Eigen::MatrixXd& Xnew) const override;
        bool _set_hyperparameters(const Eigen::VectorXd& log_lengthscale, const double& log_nugget) override;
    public:
        using Kriging::Kriging;
};

class LimitKriging : public OrdinaryKriging {
    protected:
        double _get_loocv() const override;
        void _predict(const Eigen::VectorXd& xnew, double& mean, double& sd) const override;
        Rcpp::List _predict_neumann(const Eigen::MatrixXd& Xnew) const override;
    public:
        using OrdinaryKriging::OrdinaryKriging;
};

class RationalKriging : public Kriging {
    protected:
        Eigen::VectorXd _c;
        Eigen::VectorXd _s;
        void _add_data(const Eigen::MatrixXd& Xn, const Eigen::VectorXd& yn) override;
        double _get_nllh() const override;
        double _get_loocv() const override;
        void _predict(const Eigen::VectorXd& xnew, double& mean, double& sd) const override;
        Rcpp::List _predict_neumann(const Eigen::MatrixXd& Xnew) const override;
        bool _set_hyperparameters(const Eigen::VectorXd& log_lengthscale, const double& log_nugget) override;
    public:
        RationalKriging(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, 
                        Kernel& Ker, const bool& interpolation, 
                        const std::string& loss_type);
        Eigen::VectorXd get_c() const;
};

class GeneralizedRationalKriging : public Kriging {
    protected:
        Eigen::VectorXd _c;
        Eigen::VectorXd _s;
        void _add_data(const Eigen::MatrixXd& Xn, const Eigen::VectorXd& yn) override;
        double _get_nllh() const override;
        double _get_loocv() const override;
        void _predict(const Eigen::VectorXd& xnew, double& mean, double& sd) const override;
        Rcpp::List _predict_neumann(const Eigen::MatrixXd& Xnew) const override;
        bool _set_hyperparameters(const Eigen::VectorXd& log_lengthscale, const double& log_nugget) override;
    public:
        GeneralizedRationalKriging(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, 
                                   Kernel& Ker, const bool& interpolation, 
                                   const std::string& loss_type);
        double get_c0() const;
        Eigen::VectorXd get_c() const;
};

class UniversalKriging : public Kriging {
    protected:
        std::size_t _pb;
        Rcpp::Function _bfunc;
        Eigen::VectorXd _beta;
        Eigen::MatrixXd _B;
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> _LQ;
        void _add_data(const Eigen::MatrixXd& Xn, const Eigen::VectorXd& yn) override;
        double _get_nllh() const override;
        double _get_loocv() const override;
        void _predict(const Eigen::VectorXd& xnew, double& mean, double& sd) const override;
        Rcpp::List _predict_neumann(const Eigen::MatrixXd& Xnew) const override;
        bool _set_hyperparameters(const Eigen::VectorXd& log_lengthscale, const double& log_nugget) override;
    public:
        UniversalKriging(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, 
                         Kernel& Ker, const bool& interpolation, 
                         const std::size_t& pb, Rcpp::Function bfunc, 
                         const std::string& loss_type);
        Eigen::VectorXd get_beta() const;
};

#endif
