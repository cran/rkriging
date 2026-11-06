#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include "kernel.h"
#include "kriging.h"
#include "sequential.h"

RCPP_MODULE(Kernel) {
    Rcpp::class_<Kernel>("Kernel")
        .method("get_dimension", ( std::size_t (Kernel::*)() const ) &Kernel::get_dimension)
        .method("get_lengthscale", ( Eigen::VectorXd (Kernel::*)() const ) (&Kernel::get_lengthscale))
        .method("set_lengthscale", &Kernel::set_lengthscale)
        .method("compute", ( Eigen::MatrixXd (Kernel::*)(const Eigen::MatrixXd&) const )(&Kernel::compute))
    ;

    Rcpp::class_<BaseKernel>("BaseKernel")
        .derives<Kernel>("Kernel")
    ;

    Rcpp::class_<GaussianKernel>("GaussianKernel")
        .derives<BaseKernel>("BaseKernel")
        .constructor<Eigen::VectorXd>()
        .constructor<double, Eigen::VectorXd>()
    ;

    Rcpp::class_<RQKernel>("RQKernel")
        .derives<BaseKernel>("BaseKernel")
        .constructor<Eigen::VectorXd, double>()
        .constructor<double, Eigen::VectorXd, double>()
    ;

    Rcpp::class_<IMKernel>("IMKernel")
        .derives<RQKernel>("RQKernel")
        .constructor<Eigen::VectorXd, double>()
        .constructor<double, Eigen::VectorXd, double>()
    ;

    Rcpp::class_<Matern12Kernel>("Matern12Kernel")
        .derives<BaseKernel>("BaseKernel")
        .constructor<Eigen::VectorXd>()
        .constructor<double, Eigen::VectorXd>()
    ;

    Rcpp::class_<Matern32Kernel>("Matern32Kernel")
        .derives<BaseKernel>("BaseKernel")
        .constructor<Eigen::VectorXd>()
        .constructor<double, Eigen::VectorXd>()
    ;

    Rcpp::class_<Matern52Kernel>("Matern52Kernel")
        .derives<BaseKernel>("BaseKernel")
        .constructor<Eigen::VectorXd>()
        .constructor<double, Eigen::VectorXd>()
    ;

    Rcpp::class_<MaternKernel>("MaternKernel")
        .derives<BaseKernel>("BaseKernel")
        .constructor<Eigen::VectorXd, double>()
        .constructor<double, Eigen::VectorXd, double>()
    ;

    Rcpp::class_<UDFKernel>("UDFKernel")
        .derives<BaseKernel>("BaseKernel")
        .constructor<Eigen::VectorXd, Rcpp::Function>()
        .constructor<double, Eigen::VectorXd, Rcpp::Function>()
    ;

    Rcpp::class_<MultiplicativeRQKernel>("MultiplicativeRQKernel")
        .derives<RQKernel>("RQKernel")
        .constructor<Eigen::VectorXd, double>()
        .constructor<double, Eigen::VectorXd, double>()
    ;

    Rcpp::class_<MultiplicativeIMKernel>("MultiplicativeIMKernel")
        .derives<IMKernel>("IMKernel")
        .constructor<Eigen::VectorXd, double>()
        .constructor<double, Eigen::VectorXd, double>()
    ;

    Rcpp::class_<MultiplicativeMaternKernel>("MultiplicativeMaternKernel")
        .derives<MaternKernel>("MaternKernel")
        .constructor<Eigen::VectorXd, double>()
        .constructor<double, Eigen::VectorXd, double>()
    ;

    Rcpp::class_<MultiplicativeUDFKernel>("MultiplicativeUDFKernel")
        .derives<UDFKernel>("UDFKernel")
        .constructor<Eigen::VectorXd, Rcpp::Function>()
        .constructor<double, Eigen::VectorXd, Rcpp::Function>()
    ;
}

RCPP_EXPOSED_AS(Kernel)

RCPP_MODULE(Kriging) {
    Rcpp::class_<Kriging>("Kriging")
        .method("is_interpolation", ( bool (Kriging::*)() const ) &Kriging::is_interpolation)
        .method("add_data", &Kriging::add_data)
        .method("get_data", &Kriging::get_data)
        .method("get_datasize", ( std::size_t (Kriging::*)() const ) &Kriging::get_datasize)
        .method("get_dimension", ( std::size_t (Kriging::*)() const ) &Kriging::get_dimension)
        .method("get_lengthscale", ( Eigen::VectorXd (Kriging::*)() const ) &Kriging::get_lengthscale)
        .method("get_loss_type", &Kriging::get_loss_type)
        .method("get_loss", ( double (Kriging::*)() const )(&Kriging::get_loss))
        .method("get_nllh", ( double (Kriging::*)() const )(&Kriging::get_nllh))
        .method("get_loocv", ( double (Kriging::*)() const )(&Kriging::get_loocv))
        .method("get_mu", ( double (Kriging::*)() const ) &Kriging::get_mu)
        .method("get_nu2", ( double (Kriging::*)() const ) &Kriging::get_nu2)
        .method("get_sigma2", ( double (Kriging::*)() const ) &Kriging::get_sigma2)
        .method("get_nugget", ( double (Kriging::*)() const ) &Kriging::get_nugget)
        .method("fit", ( void (Kriging::*)(bool) )(&Kriging::fit))
        .method("is_fitted", &Kriging::is_fitted)
        .method("set_lengthscale_bound", &Kriging::set_lengthscale_bound)
        .method("predict", ( Rcpp::List (Kriging::*)(const Eigen::MatrixXd&) const ) &Kriging::predict)
        .method("predict_neumann", ( Rcpp::List (Kriging::*)(const Eigen::MatrixXd&) const ) &Kriging::predict_neumann)
        .method("relax_lengthscale_constraint", &Kriging::relax_lengthscale_constraint)
        .method("set_lengthscale", &Kriging::set_lengthscale)
        .method("set_nugget", &Kriging::set_nugget)
        .method("set_nugget_bound", &Kriging::set_nugget_bound)
        .method("set_epsilon", &Kriging::set_epsilon)
        .method("set_nlopt_algorithm", &Kriging::set_nlopt_algorithm)
        .method("set_nlopt_local_algorithm", &Kriging::set_nlopt_local_algorithm)
        .method("set_nlopt_maxeval", &Kriging::set_nlopt_maxeval)
        .method("set_nlopt_population", &Kriging::set_nlopt_population)
    ;

    Rcpp::class_<OrdinaryKriging>("OrdinaryKriging")
        .derives<Kriging>("Kriging")
        .constructor<Eigen::MatrixXd, Eigen::VectorXd, Kernel&, bool, std::string>()
    ;

    Rcpp::class_<LimitKriging>("LimitKriging")
        .derives<OrdinaryKriging>("OrdinaryKriging")
        .constructor<Eigen::MatrixXd, Eigen::VectorXd, Kernel&, bool, std::string>()
    ;

    Rcpp::class_<RationalKriging>("RationalKriging")
        .derives<Kriging>("Kriging")
        .constructor<Eigen::MatrixXd, Eigen::VectorXd, Kernel&, bool, std::string>()
        .method("get_c", ( Eigen::VectorXd (RationalKriging::*)() const ) &RationalKriging::get_c)
    ;

    Rcpp::class_<GeneralizedRationalKriging>("GeneralizedRationalKriging")
        .derives<Kriging>("Kriging")
        .constructor<Eigen::MatrixXd, Eigen::VectorXd, Kernel&, bool, std::string>()
        .method("get_c0", ( double (GeneralizedRationalKriging::*)() const ) &GeneralizedRationalKriging::get_c0)
        .method("get_c", ( Eigen::VectorXd (GeneralizedRationalKriging::*)() const )
                          &GeneralizedRationalKriging::get_c)
    ;

    Rcpp::class_<UniversalKriging>("UniversalKriging")
        .derives<Kriging>("Kriging")
        .constructor<Eigen::MatrixXd, Eigen::VectorXd, Kernel&, bool, std::size_t, Rcpp::Function, std::string>()
        .method("get_beta", ( Eigen::VectorXd (UniversalKriging::*)() const ) &UniversalKriging::get_beta)
    ;
}

RCPP_EXPOSED_AS(Kriging)

RCPP_MODULE(Sequential) {
    Rcpp::class_<Sequential>("Sequential")
        .method("get_next_design", &Sequential::get_next_design)
        .method("get_next_design_from_candidates", &Sequential::get_next_design_from_candidates)
        .method("get_design_dimension", &Sequential::get_design_dimension)
        .method("get_design_lower_bound", &Sequential::get_design_lower_bound)
        .method("get_design_upper_bound", &Sequential::get_design_upper_bound)
        .method("fit_kriging", ( void (Sequential::*)(bool) )(&Sequential::fit_kriging))
        .method("add_data", &Sequential::add_data)
        .method("get_data", &Sequential::get_data)
        .method("predict", ( Rcpp::List (Sequential::*)(const Eigen::MatrixXd&) const )(&Sequential::predict))
        .method("set_design_bound", &Sequential::set_design_bound)
        .method("set_nlopt_algorithm", &Sequential::set_nlopt_algorithm)
        .method("set_nlopt_local_algorithm", &Sequential::set_nlopt_local_algorithm)
        .method("set_nlopt_maxeval", &Sequential::set_nlopt_maxeval)
        .method("set_nlopt_population", &Sequential::set_nlopt_population)
    ;

    Rcpp::class_<EmulationSequential>("EmulationSequential")
        .derives<Sequential>("Sequential")
    ;

    Rcpp::class_<OptimizationSequential>("OptimizationSequential")
        .derives<Sequential>("Sequential")
        .method("get_best_objective", &OptimizationSequential::get_best_objective)
        .method("get_best_design", &OptimizationSequential::get_best_design)
    ;

    Rcpp::class_<ALMSequential>("ALMSequential")
        .derives<EmulationSequential>("EmulationSequential")
        .constructor<Kriging&, Eigen::VectorXd, Eigen::VectorXd>()
    ;

    Rcpp::class_<UCBSequential>("UCBSequential")
        .derives<OptimizationSequential>("OptimizationSequential")
        .constructor<Kriging&, Eigen::VectorXd, Eigen::VectorXd, double>()
    ;

    Rcpp::class_<PISequential>("PISequential")
        .derives<OptimizationSequential>("OptimizationSequential")
        .constructor<Kriging&, Eigen::VectorXd, Eigen::VectorXd>()
    ;

    Rcpp::class_<EISequential>("EISequential")
        .derives<OptimizationSequential>("OptimizationSequential")
        .constructor<Kriging&, Eigen::VectorXd, Eigen::VectorXd>()
    ;

    Rcpp::class_<UDFSequential>("UDFSequential")
        .derives<Sequential>("Sequential")
        .constructor<Kriging&, Eigen::VectorXd, Eigen::VectorXd, Rcpp::Function, Rcpp::List>()
    ;
}

RCPP_EXPOSED_AS(Sequential)
