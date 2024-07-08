#ifndef KRINGING_H
#define KRINGING_H

#include <limits>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include "utils.h"
#include "kernel.h"

double nlopt_nllh(unsigned n, const double* x, double* grad, void* data);

class Kriging {
    protected:
        std::size_t n_;
        std::size_t p_;
        Eigen::MatrixXd X_;
        Eigen::VectorXd y_;
        Kernel& Ker_;
        bool interpolation_;
        double nugget_ = 1e-6;
        double epsilon_ = 1e-3;
        double mu_;
        double nu2_;
        Eigen::VectorXd a_;
        Eigen::VectorXd b_;
        Eigen::LLT<Eigen::MatrixXd> L_;
        double nllh_;
        double y_tss_;
        std::string nlopt_algorithm_ = "NLOPT_LN_SBPLX";
        std::string nlopt_local_algorithm_ = "";
        std::size_t nlopt_maxeval_ = 100;
        virtual void predict(const Eigen::VectorXd& xnew, double& mean, double& sd) = 0;
    public:
        Kriging(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, Kernel& Ker, const bool& interpolation);
        bool is_interpolation();
        virtual void add_data(const Eigen::MatrixXd& Xn, const Eigen::VectorXd& yn);
        Eigen::MatrixXd get_X();
        Eigen::VectorXd get_y();
        std::size_t get_datasize();
        std::size_t get_dimension();
        std::size_t get_lengthscale_dimension();
        double get_nllh();
        virtual double get_nllh(const Eigen::VectorXd& log_lengthscale, const double& nugget) = 0;
        void fit_hyperparameters(const Eigen::VectorXd& lengthscale_lower_bound, const Eigen::VectorXd& lengthscale_upper_bound);
        void set_lengthscale(const Eigen::VectorXd& lengthscale);
        Eigen::VectorXd get_lengthscale();
        void relax_lengthscale_constraint();
        void set_nugget(const double& nugget);
        void set_epsilon(const double& epsilon);
        void set_nlopt_algorithm(const std::string& nlopt_algorithm);
        void set_nlopt_local_algorithm(const std::string& nlopt_local_algorithm);
        void set_nlopt_maxeval(const std::size_t& nlopt_maxeval);
        virtual void set_kriging_parameters() = 0;
        void fit(const Eigen::VectorXd& lengthscale_lower_bound, const Eigen::VectorXd& lengthscale_upper_bound);
        void fit();
        double get_mu();
        double get_nu2();
        double get_sigma2();
        Eigen::VectorXd predict(const Eigen::VectorXd& xnew);
        Rcpp::List predict(const Eigen::MatrixXd& Xnew);
        Rcpp::List get_data();
        ~Kriging() {};
};

class OrdinaryKriging : public Kriging {
    protected:
        double get_nllh(const Eigen::VectorXd& log_lengthscale, const double& nugget, double& mu, double& nu2, Eigen::LLT<Eigen::MatrixXd>& L, Eigen::VectorXd& a, Eigen::VectorXd& b);
        void predict(const Eigen::VectorXd& xnew, double& mean, double& sd) override;
    public:
        using Kriging::Kriging;
        using Kriging::get_nllh;
        double get_nllh(const Eigen::VectorXd& log_lengthscale, const double& nugget) override;
        void set_kriging_parameters() override;
};

class LimitKriging : public OrdinaryKriging {
    protected:
        void predict(const Eigen::VectorXd& xnew, double& mean, double& sd) override;
    public:
        using OrdinaryKriging::OrdinaryKriging;
};

class RationalKriging : public Kriging {
    protected:
        Eigen::VectorXd c_;
        Eigen::VectorXd s_;
        double get_nllh(const Eigen::VectorXd& log_lengthscale, const double& nugget, double& mu, double& nu2, Eigen::LLT<Eigen::MatrixXd>& L, Eigen::VectorXd& a, Eigen::VectorXd& b, Eigen::VectorXd& c, Eigen::VectorXd& s);
        void predict(const Eigen::VectorXd& xnew, double& mean, double& sd) override;
    public:
        RationalKriging(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, Kernel& Ker, const bool& interpolation);
        void add_data(const Eigen::MatrixXd& Xn, const Eigen::VectorXd& yn) override;
        using Kriging::get_nllh;
        double get_nllh(const Eigen::VectorXd& log_lengthscale, const double& nugget) override;
        void set_kriging_parameters() override;
        Eigen::VectorXd get_c();
};

class GeneralizedRationalKriging : public Kriging {
    protected:
        Eigen::VectorXd c_;
        Eigen::VectorXd s_;
        double get_nllh(const Eigen::VectorXd& log_lengthscale, const double& nugget, double& mu, double& nu2, Eigen::LLT<Eigen::MatrixXd>& L, Eigen::VectorXd& a, Eigen::VectorXd& b, Eigen::VectorXd& c, Eigen::VectorXd& s);
        void predict(const Eigen::VectorXd& xnew, double& mean, double& sd) override;
    public:
        GeneralizedRationalKriging(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, Kernel& Ker, const bool& interpolation);
        void add_data(const Eigen::MatrixXd& Xn, const Eigen::VectorXd& yn) override;
        using Kriging::get_nllh;
        double get_nllh(const Eigen::VectorXd& log_lengthscale, const double& nugget) override;
        void set_kriging_parameters() override;
        double get_c0();
        Eigen::VectorXd get_c();
};

class UniversalKriging : public Kriging {
    protected:
        std::size_t pb_;
        Rcpp::Function bfunc_;
        Eigen::VectorXd beta_;
        Eigen::MatrixXd B_;
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> LQ_;
        double get_nllh(const Eigen::VectorXd& log_lengthscale, const double& nugget, double& nu2, Eigen::LLT<Eigen::MatrixXd>& L, Eigen::ColPivHouseholderQR<Eigen::MatrixXd>& LQ, Eigen::VectorXd& a, Eigen::VectorXd& b, Eigen::VectorXd& beta);
        void predict(const Eigen::VectorXd& xnew, double& mean, double& sd) override;
    public:
        UniversalKriging(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, Kernel& Ker, const bool& interpolation, const std::size_t& pb, Rcpp::Function bfunc);
        void add_data(const Eigen::MatrixXd& Xn, const Eigen::VectorXd& yn) override;
        double get_nllh(const Eigen::VectorXd& log_lengthscale, const double& nugget) override;
        void set_kriging_parameters() override;
        Eigen::VectorXd get_beta();
};

#endif
