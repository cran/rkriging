#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include "kernel.h"
#include "kriging.h"

RCPP_MODULE(Kernel) {
    Rcpp::class_<Kernel>("Kernel")
        .method("get_dimension", &Kernel::get_dimension)
        .method("get_lengthscale", ( Eigen::VectorXd (Kernel::*)() ) (&Kernel::get_lengthscale))
        .method("set_lengthscale", &Kernel::set_lengthscale)
        .method("compute", ( Eigen::MatrixXd (Kernel::*)(const Eigen::MatrixXd&) )(&Kernel::compute))
    ;

    Rcpp::class_<BaseKernel>("BaseKernel")
        .derives<Kernel>("Kernel")
    ;
    
    Rcpp::class_<GaussianKernel>("GaussianKernel")
        .derives<BaseKernel>("BaseKernel")
        .constructor<Eigen::VectorXd>()
        .constructor<double,Eigen::VectorXd>()
    ;

    Rcpp::class_<RQKernel>("RQKernel")
        .derives<BaseKernel>("BaseKernel")
        .constructor<Eigen::VectorXd,double>()
        .constructor<double,Eigen::VectorXd,double>()
    ;

    Rcpp::class_<Matern12Kernel>("Matern12Kernel")
        .derives<BaseKernel>("BaseKernel")
        .constructor<Eigen::VectorXd>()
        .constructor<double,Eigen::VectorXd>()
    ;

    Rcpp::class_<Matern32Kernel>("Matern32Kernel")
        .derives<BaseKernel>("BaseKernel")
        .constructor<Eigen::VectorXd>()
        .constructor<double,Eigen::VectorXd>()
    ;

    Rcpp::class_<Matern52Kernel>("Matern52Kernel")
        .derives<BaseKernel>("BaseKernel")
        .constructor<Eigen::VectorXd>()
        .constructor<double,Eigen::VectorXd>()
    ;

    Rcpp::class_<MaternKernel>("MaternKernel")
        .derives<BaseKernel>("BaseKernel")
        .constructor<Eigen::VectorXd,double>()
        .constructor<double,Eigen::VectorXd,double>()
    ;

    Rcpp::class_<UDFKernel>("UDFKernel")
        .derives<BaseKernel>("BaseKernel")
        .constructor<Eigen::VectorXd,Rcpp::Function>()
        .constructor<double,Eigen::VectorXd,Rcpp::Function>()
    ;

    Rcpp::class_<MultiplicativeRQKernel>("MultiplicativeRQKernel")
        .derives<RQKernel>("RQKernel")
        .constructor<Eigen::VectorXd,double>()
        .constructor<double,Eigen::VectorXd,double>()
    ;

    Rcpp::class_<MultiplicativeMaternKernel>("MultiplicativeMaternKernel")
        .derives<MaternKernel>("MaternKernel")
        .constructor<Eigen::VectorXd,double>()
        .constructor<double,Eigen::VectorXd,double>()
    ;

    Rcpp::class_<MultiplicativeUDFKernel>("MultiplicativeUDFKernel")
        .derives<UDFKernel>("UDFKernel")
        .constructor<Eigen::VectorXd,Rcpp::Function>()
        .constructor<double,Eigen::VectorXd,Rcpp::Function>()
    ;
}

RCPP_EXPOSED_AS(Kernel)

RCPP_MODULE(Kriging) {
    Rcpp::class_<Kriging>("Kriging")
        .method("is_interpolation", &Kriging::is_interpolation)
        .method("get_nllh", ( double (Kriging::*)() )(&Kriging::get_nllh))
        .method("add_data", &Kriging::add_data)
        .method("get_data", &Kriging::get_data)
        .method("get_datasize", &Kriging::get_datasize)
        .method("get_dimension", &Kriging::get_dimension)
        .method("get_lengthscale", &Kriging::get_lengthscale)
        .method("set_lengthscale", &Kriging::set_lengthscale)
        .method("relax_lengthscale_constraint", &Kriging::relax_lengthscale_constraint)
        .method("set_nugget", &Kriging::set_nugget)
        .method("set_epsilon", &Kriging::set_epsilon)
        .method("set_nlopt_algorithm", &Kriging::set_nlopt_algorithm)
        .method("set_nlopt_local_algorithm", &Kriging::set_nlopt_local_algorithm)
        .method("set_nlopt_maxeval", &Kriging::set_nlopt_maxeval)
        .method("set_kriging_parameters", &Kriging::set_kriging_parameters)
        .method("fit", ( void (Kriging::*)(const Eigen::VectorXd&,const Eigen::VectorXd&) )(&Kriging::fit))
        .method("fit", ( void (Kriging::*)() )(&Kriging::fit))
        .method("get_mu", &Kriging::get_mu)
        .method("get_nu2", &Kriging::get_nu2)
        .method("get_sigma2", &Kriging::get_sigma2)
        .method("predict", ( Rcpp::List (Kriging::*)(const Eigen::MatrixXd&) ) &Kriging::predict)
    ;

    Rcpp::class_<OrdinaryKriging>("OrdinaryKriging")
        .derives<Kriging>("Kriging")
        .constructor<Eigen::MatrixXd,Eigen::VectorXd,Kernel&,bool>()
    ;

    Rcpp::class_<LimitKriging>("LimitKriging")
        .derives<OrdinaryKriging>("OrdinaryKriging")
        .constructor<Eigen::MatrixXd,Eigen::VectorXd,Kernel&,bool>()
    ;

    Rcpp::class_<RationalKriging>("RationalKriging")
        .derives<Kriging>("Kriging")
        .constructor<Eigen::MatrixXd,Eigen::VectorXd,Kernel&,bool>()
        .method("get_c", &RationalKriging::get_c)
    ;

    Rcpp::class_<GeneralizedRationalKriging>("GeneralizedRationalKriging")
        .derives<Kriging>("Kriging")
        .constructor<Eigen::MatrixXd,Eigen::VectorXd,Kernel&,bool>()
        .method("get_c0", &GeneralizedRationalKriging::get_c0)
        .method("get_c", &GeneralizedRationalKriging::get_c)
    ;

    Rcpp::class_<UniversalKriging>("UniversalKriging")
        .derives<Kriging>("Kriging")
        .constructor<Eigen::MatrixXd,Eigen::VectorXd,Kernel&,bool,std::size_t,Rcpp::Function>()
        .method("get_beta", &UniversalKriging::get_beta)
    ;
}

RCPP_EXPOSED_AS(Kriging)
