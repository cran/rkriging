#' Kriging
#'
#' @useDynLib Kriging, .registration = TRUE
#' @import Rcpp, RcppEigen, BH, nloptr
#' 
#' @noRd
Rcpp::loadModule("Kriging", TRUE)

#' @noRd
Get.Lengthscale <- function(X, lengthscale=list(value=NULL,lower.bound=NULL,upper.bound=NULL)) {
  n <- nrow(X)
  p <- ncol(X)
  # check input initial / lower bound / upper bound
  if (!all(sapply(lengthscale, function(x) is.null(x)|length(x)==p)))
    stop(sprintf("Kernel lengthscale (initial/lower bound/upper bound) dimension must match with the problem dimension (%d).", p))
  if (any(lengthscale[["lower.bound"]]>lengthscale[["upper.bound"]]))
    stop("Kernel lengthscale lower bound can not be greater than its upper bound.")
  if (any(lengthscale[["value"]]<lengthscale[["lower.bound"]])|any(lengthscale[["value"]]>lengthscale[["upper.bound"]]))
    stop("Kernel lengthscale initial must be between lower bound and upper bound.")
  # compute initial / lower bound / upper bound if needed
  if (any(sapply(lengthscale, function(x) is.null(x)))) {
    X.pdist <- dist(X)
    X.pdist.med <- median(X.pdist)
    X.pdist <- as.matrix(X.pdist)
    X.ndist.max <- max(sapply(1:n, function(i) min(X.pdist[i,-i]))) 
    ls.ini <- rep(max(X.pdist.med, 1.1*X.ndist.max/2), p) / sqrt(2)
    ls.lb <- rep(X.ndist.max/2, p) / sqrt(2)
    ls.ub <- pmax(ls.ini*(ls.ini / ls.lb), 100)
    if (is.null(lengthscale[["value"]])) {
      lengthscale[["value"]] <- ls.ini
      if (!is.null(lengthscale[["lower.bound"]]) & !is.null(lengthscale[["upper.bound"]])) {
        for (i in 1:p)
          if (lengthscale[["value"]][i]<lengthscale[["lower.bound"]][i] | lengthscale[["value"]][i]>lengthscale[["upper.bound"]][i])
            lengthscale[["value"]][i] <- 0.5 * (lengthscale[["lower.bound"]][i] + lengthscale[["upper.bound"]][i])
      } else {
        if (!is.null(lengthscale[["lower.bound"]]))
          lengthscale[["value"]] <- pmax(lengthscale[["value"]], lengthscale[["lower.bound"]])
        if (!is.null(lengthscale[["upper.bound"]]))
          lengthscale[["value"]] <- pmin(lengthscale[["value"]], lengthscale[["upper.bound"]])
      }
    }
    if (is.null(lengthscale[["lower.bound"]]))
      lengthscale[["lower.bound"]] <- pmin(lengthscale[["value"]], ls.lb)
    if (is.null(lengthscale[["upper.bound"]]))
      lengthscale[["upper.bound"]] <- pmax(lengthscale[["value"]], ls.ub)
  }
  return (lengthscale)
}

#' @noRd
Set.Kriging.NLOPT.Parameters <- function(kriging, nlopt.parameters=list()) {
  if (!grepl("^Rcpp.*Kriging$", class(kriging)))
    stop("An invalid kriging object is provided.")
  # rkriging supported optimization algorithm available from nlopt
  global.optimization.algorithms <- c(
    "NLOPT_GN_DIRECT", "NLOPT_GN_DIRECT_L", "NLOPT_GN_DIRECT_L_RAND", 
    "NLOPT_GN_DIRECT_NOSCAL","NLOPT_GN_DIRECT_L_NOSCAL", "NLOPT_GN_DIRECT_L_RAND_NOSCAL", 
    "NLOPT_GN_ORIG_DIRECT", "NLOPT_GN_ORIG_DIRECT_L", "NLOPT_GN_ISRES",
    "NLOPT_GN_CRS2_LM", "NLOPT_GN_MLSL", "NLOPT_GN_MLSL_LDS"
  )
  local.optimization.algorithms <- c(
    "NLOPT_LN_COBYLA", "NLOPT_LN_NEWUOA", "NLOPT_LN_BOBYQA",
    "NLOPT_LN_NELDERMEAD", "NLOPT_LN_SBPLX", "NLOPT_LN_AUGLAG"
  )
  optimization.algorithms <- c(global.optimization.algorithms, local.optimization.algorithms)
  # get nlopt parameters 
  if ("algorithm" %in% names(nlopt.parameters)) {
    nlopt.algorithm <- nlopt.parameters[["algorithm"]]
    if (!(nlopt.algorithm %in% optimization.algorithms)) 
      stop(sprintf("%s optimization algorithm is not supported. Only gradient-free optimization is supported: %s", 
                   nlopt.algorithm, paste(optimization.algorithms, collapse=", ")))
    kriging$set_nlopt_algorithm(nlopt.algorithm)
    if (nlopt.algorithm %in% c("NLOPT_GN_MLSL","NLOPT_GN_MLSL_LDS")) {
      if ("local.algorithm" %in% names(nlopt.parameters)) {
        nlopt.local.algorithm <- nlopt.parameters[["local.algorithm"]]
        if (!(nlopt.local.algorithm %in% local.optimization.algorithms))
          stop(sprintf("%s optimization algorithm is not supported. Only local gradient-free optimization is supported: %s", 
                       nlopt.local.algorithm, paste(local.optimization.algorithms, collapse=", ")))
        kriging$set_nlopt_local_algorithm(nlopt.local.algorithm)
      } else {
        kriging$set_nlopt_local_algorithm("NLOPT_LN_SBPLX")
      }
    }
  } 
  if ("maxeval" %in% names(nlopt.parameters)) {
    nlopt.maxeval <- as.integer(nlopt.parameters[["maxeval"]])
    if (nlopt.maxeval <= 0) 
      stop(sprintf("nlopt optimization maximum evaluation must be a positive integer."))
    kriging$set_nlopt_maxeval(nlopt.maxeval)
  }
  return (kriging)
}

#' @title 
#' Fit or Create Kriging Models
#' 
#' @description
#' This function provides a common interface to fit various kriging models from the data.
#' OK (\link{Ordinary.Kriging}), UK (\link{Universal.Kriging}), LK (\link{Limit.Kriging}), 
#' RK (\link{Rational.Kriging}), and GRK (\link{Generalized.Rational.Kriging}) are supported in this package.
#' 
#' @details 
#' Kriging gives the best linear unbiased predictor given the data. 
#' It can also be given a Bayesian interpretation by assuming a Gaussian process (GP) prior for the underlying function.
#' Please see Santner et al. (2003), Rasmussen and Williams (2006), and Joseph (2024) for details.
#' 
#' For data from deterministic computer experiments, use \code{interpolation=TRUE} and will give an interpolator. 
#' For noisy data, use \code{interpolation=FALSE}, which will give an approximator of the underlying function.
#' Currently, approximation is supported for OK (\link{Ordinary.Kriging}) and UK (\link{Universal.Kriging}). 
#' 
#' The kernel choices are required and can be specified by 
#' (i) providing the kernel class object to \code{kernel}
#' or (ii) specifying the kernel type and other parameters in \code{kernel.parameters}. 
#' Please see examples section for detail usages. 
#' 
#' When the lengthscale / correlation parameters are known, simply provide them in 
#' (i) the kernel class object to \code{kernel} or (ii) in \code{kernel.parameters},
#' and set \code{fit=FALSE}. The kriging model will be fitted with the user provided parameters.
#'  
#' When the lengthscale / correlation parameters are unknown, 
#' they can be estimated via Maximum Likelihood method by setting \code{fit=TRUE}. 
#' The initial / lower bound / upper bound of the lengthscale parameters can be provided in \code{kernel.parameters}, 
#' otherwise a good initial and range would be estimated from the data. 
#' The optimization is performed via \href{https://nlopt.readthedocs.io/en/latest/}{NLopt}, 
#' a open-source library for nonlinear optimization. 
#' All gradient-free optimization methods in \href{https://nlopt.readthedocs.io/en/latest/}{NLopt} 
#' are supported and can be specified in \code{nlopt.parameters}.
#' See \code{nloptr::nloptr.print.options()} for the list of available derivative-free algorithms (prefix with NLOPT_GN or NLOPT_LN).  
#' The maximum number of optimization steps can also be defined in \code{nlopt.parameters}.
#' Please see examples section for detail usages. 
#' 
#' @param X a matrix for input (feature)
#' @param y a vector for output (target), only one-dimensional output is supported
#' @param interpolation whether to interpolate, for noisy data please set \code{interpolate=FALSE}
#' @param fit whether to fit the length scale parameters from data 
#' @param model choice of kriging model: OK, UK, LK, RK, GRK
#' @param model.parameters a list of parameters required for the specific kriging model, e.g. basis functions for universal kriging
#' @param kernel a kernel class object
#' @param kernel.parameters a list of parameters required for the kernel, if no kernel class object is provided
#' @param nlopt.parameters a list of parameters required for NLopt, including choice of optimization algorithm and maximum number of evaluation
#' 
#' @return 
#' A Kriging Class Object.
#' 
#' @author 
#' Chaofan Huang and V. Roshan Joseph
#' 
#' @references
#' Joseph, V. R. (2006). \emph{Limit kriging}. Technometrics, 48(4), 458-466.
#' 
#' Joseph, V. R. (2024). Rational Kriging. \emph{Journal of the American Statistical Association}.
#' 
#' Rasmussen, C. E. & Williams, C. K. (2006). \emph{Gaussian Processes for Machine Learning}. The MIT Press.
#' 
#' Santner, T. J., Williams, B. J., Notz, W. I., & Williams, B. J. (2003). \emph{The design and analysis of computer experiments (Vol. 1)}. New York: Springer.
#' 
#' @seealso 
#' \link{Predict.Kriging}, \link{Get.Kriging.Parameters}, \link{Get.Kernel}, 
#' \link{Ordinary.Kriging}, \link{Universal.Kriging}, \link{Limit.Kriging}, 
#' \link{Rational.Kriging}, \link{Generalized.Rational.Kriging}.
#' 
#' @export
#' 
#' @examples
#' # one dimensional example 
#' f <- function(x) {
#'   x <- 0.5 + 2*x
#'   y <- sin(10*pi*x)/(2*x) + (x-1)^4
#'   return (y)
#' }
#' 
#' set.seed(1234)
#' # train set
#' n <- 30
#' p <- 1
#' X <- matrix(runif(n),ncol=p)
#' y <- apply(X, 1, f)
#  # test set
#' newX <- matrix(seq(0,1,length=1001), ncol=p)
#' 
#' #############################################################################
#' ################ Minimal Example for Fitting a Kriging Model ################
#' #############################################################################
#' # Ordinary Kriging
#' kriging <- Fit.Kriging(X, y, interpolation=TRUE, fit=TRUE, model="OK",
#'                        kernel.parameters=list(type="Gaussian"))
#' pred <- Predict.Kriging(kriging, newX)
#' plot(newX, f(newX), "l")
#' points(X, y, pch=16, col="blue")
#' lines(newX, pred$mean, col="red", lty=2)
#' lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#' lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#' Get.Kriging.Parameters(kriging)
#' 
#' # Universal Kriging
#' basis.function <- function(x) {c(1,x[1],x[1]^2)}
#' kriging <- Fit.Kriging(X, y, interpolation=TRUE, fit=TRUE, model="UK",
#'                        model.parameters=list(basis.function=basis.function),
#'                        kernel.parameters=list(type="Gaussian"))
#' pred <- Predict.Kriging(kriging, newX)
#' plot(newX, f(newX), "l")
#' points(X, y, pch=16, col="blue")
#' lines(newX, pred$mean, col="red", lty=2)
#' lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#' lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#' Get.Kriging.Parameters(kriging)
#' 
#' # Limit Kriging
#' kriging <- Fit.Kriging(X, y, interpolation=TRUE, fit=TRUE, model="LK",
#'                        kernel.parameters=list(type="Gaussian"))
#' pred <- Predict.Kriging(kriging, newX)
#' plot(newX, f(newX), "l")
#' points(X, y, pch=16, col="blue")
#' lines(newX, pred$mean, col="red", lty=2)
#' lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#' lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#' Get.Kriging.Parameters(kriging)
#' 
#' # Rational Kriging
#' kriging <- Fit.Kriging(X, y, interpolation=TRUE, fit=TRUE, model="RK",
#'                        kernel.parameters=list(type="RQ",alpha=1))
#' pred <- Predict.Kriging(kriging, newX)
#' plot(newX, f(newX), "l")
#' points(X, y, pch=16, col="blue")
#' lines(newX, pred$mean, col="red", lty=2)
#' lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#' lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#' Get.Kriging.Parameters(kriging)
#' 
#' # Generalized Rational Kriging
#' kriging <- Fit.Kriging(X, y, interpolation=TRUE, fit=TRUE, model="GRK",
#'                        kernel.parameters=list(type="RQ",alpha=1))
#' pred <- Predict.Kriging(kriging, newX)
#' plot(newX, f(newX), "l")
#' points(X, y, pch=16, col="blue")
#' lines(newX, pred$mean, col="red", lty=2)
#' lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#' lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#' Get.Kriging.Parameters(kriging)
#' 
#' #############################################################################
#' ################ Fitting a Kriging Model with Kernel Object #################
#' #############################################################################
#' kernel <- Gaussian.Kernel(rep(1,p))
#' kriging <- Fit.Kriging(X, y, interpolation=TRUE, fit=TRUE, model="OK", kernel=kernel)
#' pred <- Predict.Kriging(kriging, newX)
#' plot(newX, f(newX), "l")
#' points(X, y, pch=16, col="blue")
#' lines(newX, pred$mean, col="red", lty=2)
#' lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#' lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#' Get.Kriging.Parameters(kriging)
#' 
#' #############################################################################
#' ############### Creating a Kriging Model with Kernel Object #################
#' #############################################################################
#' # set fit = FALSE to create Kriging model with user provided kernel
#' # no optimization for the length scale parameters 
#' kernel <- Gaussian.Kernel(rep(1e-1,p))
#' kriging <- Fit.Kriging(X, y, interpolation=TRUE, fit=FALSE, model="OK", kernel=kernel)
#' pred <- Predict.Kriging(kriging, newX)
#' plot(newX, f(newX), "l")
#' points(X, y, pch=16, col="blue")
#' lines(newX, pred$mean, col="red", lty=2)
#' lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#' lines(newX, pred$mean+2*pred$sd, col="red", lty=3) 
#' Get.Kriging.Parameters(kriging)
#' 
#' #############################################################################
#' ############ Fitting a Kriging Model with Range for Lengthscale #############
#' #############################################################################
#' kernel.parameters <- list(
#'     type = 'Gaussian',
#'     lengthscale = rep(1,p), # initial value
#'     lengthscale.lower.bound = rep(1e-3,p), # lower bound
#'     lengthscale.upper.bound = rep(1e3,p) # upper bound
#' ) # if not provided, a good estimate would be computed from data
#' kriging <- Fit.Kriging(X, y, interpolation=TRUE, fit=TRUE, model="OK",
#'                        kernel.parameters=kernel.parameters)
#' pred <- Predict.Kriging(kriging, newX)
#' plot(newX, f(newX), "l")
#' points(X, y, pch=16, col="blue")
#' lines(newX, pred$mean, col="red", lty=2)
#' lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#' lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#' Get.Kriging.Parameters(kriging)
#' 
#' #############################################################################
#' ########### Fitting a Kriging Model with Different Optimization #############
#' #############################################################################
#' nlopt.parameters <- list(
#'     algorithm = 'NLOPT_GN_DIRECT', # optimization method
#'     maxeval = 250 # maximum number of evaluation
#' )
#' kriging <- Fit.Kriging(X, y, interpolation=TRUE, fit=TRUE, model="OK",
#'                        kernel.parameters=list(type="Gaussian"),
#'                        nlopt.parameters=nlopt.parameters)
#' pred <- Predict.Kriging(kriging, newX)
#' plot(newX, f(newX), "l")
#' points(X, y, pch=16, col="blue")
#' lines(newX, pred$mean, col="red", lty=2)
#' lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#' lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#' Get.Kriging.Parameters(kriging)
#' 
#' # Global-Local Optimization
#' nlopt.parameters <- list(
#'     algorithm = 'NLOPT_GN_MLSL_LDS', # optimization method
#'     local.algorithm = 'NLOPT_LN_SBPLX', # local algorithm
#'     maxeval = 250 # maximum number of evaluation
#' )
#' kriging <- Fit.Kriging(X, y, interpolation=TRUE, fit=TRUE, model="OK",
#'                        kernel.parameters=list(type="Gaussian"),
#'                        nlopt.parameters=nlopt.parameters)
#' pred <- Predict.Kriging(kriging, newX)
#' plot(newX, f(newX), "l")
#' points(X, y, pch=16, col="blue")
#' lines(newX, pred$mean, col="red", lty=2)
#' lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#' lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#' Get.Kriging.Parameters(kriging)
#' 
#' #############################################################################
#' ################# Fitting a Kriging Model from Noisy Data ###################
#' #############################################################################
#' y <- y + 0.1 * rnorm(length(y))
#' kriging <- Fit.Kriging(X, y, interpolation=FALSE, fit=TRUE, model="OK",
#'                        kernel.parameters=list(type="Gaussian"))
#' pred <- Predict.Kriging(kriging, newX)
#' plot(newX, f(newX), "l")
#' points(X, y, pch=16, col="blue")
#' lines(newX, pred$mean, col="red", lty=2)
#' lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#' lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#' Get.Kriging.Parameters(kriging)
#' 
Fit.Kriging <- function(X, y, interpolation=TRUE, fit=TRUE,
                        model="OK", model.parameters=list(), 
                        kernel=NULL, kernel.parameters=list(), 
                        nlopt.parameters=list()) {
  # check inputs
  if (!is.matrix(X)) X <- as.matrix(X)
  if (is.matrix(y)) y <- c(y)
  if (!is.character(model)) 
    stop("kriging model must be a string.")
  # get lengthscale parameters if not provided
  lengthscale <- list(
    value = kernel.parameters[["lengthscale"]],
    lower.bound = kernel.parameters[["lengthscale.lower.bound"]],
    upper.bound = kernel.parameters[["lengthscale.upper.bound"]]
  )
  if (fit) lengthscale <- Get.Lengthscale(X, lengthscale)
  # get kernel class
  if (is.null(kernel)) {
    if (!("type" %in% names(kernel.parameters))) 
      stop("either kernel object or type of kernel must be provided in kernel.parameters.")
    if (!("value" %in% names(lengthscale)))
      stop("either kernel object or lengthscale of kernel must be provided in kernel.parameters when fit = FALSE.")
    kernel <- Get.Kernel(lengthscale$value, kernel.parameters[["type"]], kernel.parameters)
  } else {
    if (!grepl("^Rcpp.*Kernel$", class(kernel)))
      stop("An invalid kernel object is provided.")
    if (ncol(X) != kernel$get_dimension())
      stop(sprintf("Dimension of X (%d) does not match with the kernel dimension (%d).", 
                   ncol(X), kernel$get_dimension()))
  }
  # check if noisy data is supported
  if (!interpolation & model %in% c("LK","RK","GRK")) {
    warning(sprintf("Only interpolation is supported for %s. Setting interpolation to TRUE.", model))
    interpolation <- TRUE
  }
  # initialize kriging class
  if (model == "OK") {
    kriging <- new(OrdinaryKriging, X, y, kernel, interpolation)
  } else if (model == "LK") {
    kriging <- new(LimitKriging, X, y, kernel, interpolation)
  } else if (model == "RK") {
    kriging <- new(RationalKriging, X, y, kernel, interpolation)
  } else if (model == "GRK") {
    kriging <- new(GeneralizedRationalKriging, X, y, kernel, interpolation)
  } else if (model == "UK") {
    if (!("basis.function" %in% names(model.parameters)))
      stop("basis.function must be provided for Universal Kriging (UK).")
    basis.function <- model.parameters[["basis.function"]]
    if (!is.function(basis.function))
      stop("basis.function must a function for Universal Kriging (UK).")
    no.basis.function <- length(basis.function(X[1,]))
    kriging <- new(UniversalKriging, X, y, kernel, interpolation, no.basis.function, basis.function)
  } else {
    stop(sprintf("No Kriging class %s.", model))
  }
  if (fit) {
    # set nlopt parameters 
    kriging <- Set.Kriging.NLOPT.Parameters(kriging, nlopt.parameters)
    # fit kriging
    kriging$fit(lengthscale$lower.bound, lengthscale$upper.bound)
  } else {
    # compute kriging parameters for the given lengthscale
    kriging$set_kriging_parameters()
  }
  return (kriging)
}

#' @title 
#' Ordinary Kriging
#' 
#' @description
#' This functions fits the ordinary kriging model to the data.
#' 
#' @details
#' Ordinary kriging assumes a constant mean. Please see Santner et al. (2003) for details.
#' 
#' For data from deterministic computer experiments, use \code{interpolation=TRUE} and will give an interpolator. 
#' For noisy data, use \code{interpolation=FALSE}, which will give an approximator of the underlying function.
#' 
#' The kernel choices are required and can be specified by 
#' (i) providing the kernel class object to \code{kernel}
#' or (ii) specifying the kernel type and other parameters in \code{kernel.parameters}. 
#' Please see examples section of \link{Fit.Kriging} for detail usages. 
#' 
#' When the lengthscale / correlation parameters are unknown, 
#' all parameters including the constant mean can be estimated via Maximum Likelihood method by setting \code{fit=TRUE}. 
#' The initial / lower bound / upper bound of the lengthscale parameters can be provided in \code{kernel.parameters}, 
#' otherwise a good initial and range would be estimated from the data. 
#' The optimization is performed via \href{https://nlopt.readthedocs.io/en/latest/}{NLopt}, 
#' a open-source library for nonlinear optimization. 
#' All gradient-free optimization methods in \href{https://nlopt.readthedocs.io/en/latest/}{NLopt} 
#' are supported and can be specified in \code{nlopt.parameters}.
#' See \code{nloptr::nloptr.print.options()} for the list of available derivative-free algorithms (prefix with NLOPT_GN or NLOPT_LN).  
#' The maximum number of optimization steps can also be defined in \code{nlopt.parameters}.
#' Please see examples section of \link{Fit.Kriging} for detail usages.
#' 
#' @param X a matrix for input (feature)
#' @param y a vector for output (target), only one-dimensional output is supported
#' @param interpolation interpolation whether to interpolate, for noisy data please set \code{interpolate=FALSE}
#' @param fit whether to fit the length scale parameters from data 
#' @param kernel a kernel class object 
#' @param kernel.parameters a list of parameters required for the kernel, if no kernel class object is provided
#' @param nlopt.parameters a list of parameters required for NLopt, including choice of optimization algorithm and maximum number of evaluation
#' 
#' @return 
#' A Ordinary Kriging Class Object.
#' 
#' @author 
#' Chaofan Huang and V. Roshan Joseph
#' 
#' @references
#' Santner, T. J., Williams, B. J., Notz, W. I., & Williams, B. J. (2003). \emph{The design and analysis of computer experiments (Vol. 1)}. New York: Springer.
#' 
#' @seealso 
#' \link{Fit.Kriging}, \link{Predict.Kriging}, \link{Get.Kriging.Parameters}.
#' 
#' @export
#' 
#' @examples
#' # one dimensional example 
#' f <- function(x) {
#'   x <- 0.5 + 2*x
#'   y <- sin(10*pi*x)/(2*x) + (x-1)^4
#'   return (y)
#' }
#' 
#' set.seed(1234)
#' # train set
#' n <- 30
#' p <- 1
#' X <- matrix(runif(n),ncol=p)
#' y <- apply(X, 1, f)
#  # test set
#' newX <- matrix(seq(0,1,length=1001), ncol=p)
#' 
#' # approach 1
#' kriging <- Ordinary.Kriging(X, y, interpolation=TRUE, fit=TRUE, 
#'                             kernel.parameters=list(type="Gaussian"))
#' pred <- Predict.Kriging(kriging, newX)
#' plot(newX, f(newX), "l")
#' points(X, y, pch=16, col="blue")
#' lines(newX, pred$mean, col="red", lty=2)
#' lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#' lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#' Get.Kriging.Parameters(kriging)
#' 
#' # approach 2
#' kriging <- Fit.Kriging(X, y, interpolation=TRUE, fit=TRUE, model="OK",
#'                        kernel.parameters=list(type="Gaussian"))
#' pred <- Predict.Kriging(kriging, newX)
#' plot(newX, f(newX), "l")
#' points(X, y, pch=16, col="blue")
#' lines(newX, pred$mean, col="red", lty=2)
#' lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#' lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#' Get.Kriging.Parameters(kriging)
#' 
Ordinary.Kriging <- function(X, y, interpolation=TRUE, fit=TRUE,
                            kernel=NULL, kernel.parameters=list(), 
                            nlopt.parameters=list()) {
  # check inputs
  if (!is.matrix(X)) X <- as.matrix(X)
  if (is.matrix(y)) y <- c(y)
  # get lengthscale parameters if not provided
  lengthscale <- list(
    value = kernel.parameters[["lengthscale"]],
    lower.bound = kernel.parameters[["lengthscale.lower.bound"]],
    upper.bound = kernel.parameters[["lengthscale.upper.bound"]]
  )
  if (fit) lengthscale <- Get.Lengthscale(X, lengthscale)
  # get kernel class
  if (is.null(kernel)) {
    if (!("type" %in% names(kernel.parameters))) 
      stop("either kernel object or type of kernel must be provided in kernel.parameters.")
    if (!("value" %in% names(lengthscale)))
      stop("either kernel object or lengthscale of kernel must be provided in kernel.parameters when fit = FALSE.")
    kernel <- Get.Kernel(lengthscale$value, kernel.parameters[["type"]], kernel.parameters)
  } else {
    if (!grepl("^Rcpp.*Kernel$", class(kernel)))
      stop("An invalid kernel object is provided.")
    if (ncol(X) != kernel$get_dimension())
      stop(sprintf("Dimension of X (%d) does not match with the kernel dimension (%d).", 
                   ncol(X), kernel$get_dimension()))
  }
  # initialize kriging class
  kriging <- new(OrdinaryKriging, X, y, kernel, interpolation)
  if (fit) {
    # set nlopt parameters 
    kriging <- Set.Kriging.NLOPT.Parameters(kriging, nlopt.parameters)
    # fit kriging
    kriging$fit(lengthscale$lower.bound, lengthscale$upper.bound)
  } else {
    # compute kriging parameters for the given lengthscale
    kriging$set_kriging_parameters()
  }
  return (kriging)
}

#' @title 
#' Limit Kriging
#' 
#' @description
#' This functions fits the limit kriging model to the data.
#' 
#' @details 
#' Limit kriging avoids the mean reversion issue of ordinary kriging. Please see Joseph (2006) for details.
#' Only interpolation is available. Noisy output is not supported for limit kriging.  
#' 
#' The kernel choices are required and can be specified by 
#' (i) providing the kernel class object to \code{kernel}
#' or (ii) specifying the kernel type and other parameters in \code{kernel.parameters}. 
#' Please see examples section of \link{Fit.Kriging} for detail usages. 
#' 
#' When the lengthscale / correlation parameters are unknown, 
#' all parameters including the constant mean can be estimated via Maximum Likelihood method by setting \code{fit=TRUE}. 
#' The initial / lower bound / upper bound of the lengthscale parameters can be provided in \code{kernel.parameters}, 
#' otherwise a good initial and range would be estimated from the data. 
#' The optimization is performed via \href{https://nlopt.readthedocs.io/en/latest/}{NLopt}, 
#' a open-source library for nonlinear optimization. 
#' All gradient-free optimization methods in \href{https://nlopt.readthedocs.io/en/latest/}{NLopt} 
#' are supported and can be specified in \code{nlopt.parameters}.
#' See \code{nloptr::nloptr.print.options()} for the list of available derivative-free algorithms (prefix with NLOPT_GN or NLOPT_LN).  
#' The maximum number of optimization steps can also be defined in \code{nlopt.parameters}.
#' Please see examples section of \link{Fit.Kriging} for detail usages.
#' 
#' @param X a matrix for input (feature)
#' @param y a vector for output (target), only one-dimensional output is supported
#' @param fit whether to fit the length scale parameters from data 
#' @param kernel a kernel class object 
#' @param kernel.parameters a list of parameters required for the kernel, if no kernel class object is provided
#' @param nlopt.parameters a list of parameters required for NLopt, including choice of optimization algorithm and maximum number of evaluation
#' 
#' @return 
#' A Limit Kriging Class Object.
#' 
#' @author 
#' Chaofan Huang and V. Roshan Joseph
#' 
#' @references
#' Joseph, V. R. (2006). \emph{Limit kriging}. Technometrics, 48(4), 458-466.
#' 
#' @seealso 
#' \link{Fit.Kriging}, \link{Predict.Kriging}, \link{Get.Kriging.Parameters}.
#' 
#' @export
#' 
#' @examples
#' # one dimensional example 
#' f <- function(x) {
#'   x <- 0.5 + 2*x
#'   y <- sin(10*pi*x)/(2*x) + (x-1)^4
#'   return (y)
#' }
#' 
#' set.seed(1234)
#' # train set
#' n <- 30
#' p <- 1
#' X <- matrix(runif(n),ncol=p)
#' y <- apply(X, 1, f)
#  # test set
#' newX <- matrix(seq(0,1,length=1001), ncol=p)
#' 
#' # approach 1
#' kriging <- Limit.Kriging(X, y, fit=TRUE, kernel.parameters=list(type="Gaussian"))
#' pred <- Predict.Kriging(kriging, newX)
#' plot(newX, f(newX), "l")
#' points(X, y, pch=16, col="blue")
#' lines(newX, pred$mean, col="red", lty=2)
#' lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#' lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#' Get.Kriging.Parameters(kriging)
#' 
#' # approach 2
#' kriging <- Fit.Kriging(X, y, interpolation=TRUE, fit=TRUE, model="LK",
#'                        kernel.parameters=list(type="Gaussian"))
#' pred <- Predict.Kriging(kriging, newX)
#' plot(newX, f(newX), "l")
#' points(X, y, pch=16, col="blue")
#' lines(newX, pred$mean, col="red", lty=2)
#' lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#' lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#' Get.Kriging.Parameters(kriging)
#' 
Limit.Kriging <- function(X, y, fit=TRUE, 
                          kernel=NULL, kernel.parameters=list(), 
                          nlopt.parameters=list()) {
  # only interpolation is supported for the first release
  interpolation <- TRUE
  # check inputs
  if (!is.matrix(X)) X <- as.matrix(X)
  if (is.matrix(y)) y <- c(y)
  # get lengthscale parameters if not provided
  lengthscale <- list(
    value = kernel.parameters[["lengthscale"]],
    lower.bound = kernel.parameters[["lengthscale.lower.bound"]],
    upper.bound = kernel.parameters[["lengthscale.upper.bound"]]
  )
  if (fit) lengthscale <- Get.Lengthscale(X, lengthscale)
  # get kernel class
  if (is.null(kernel)) {
    if (!("type" %in% names(kernel.parameters))) 
      stop("either kernel object or type of kernel must be provided in kernel.parameters.")
    if (!("value" %in% names(lengthscale)))
      stop("either kernel object or lengthscale of kernel must be provided in kernel.parameters when fit = FALSE.")
    kernel <- Get.Kernel(lengthscale$value, kernel.parameters[["type"]], kernel.parameters)
  } else {
    if (!grepl("^Rcpp.*Kernel$", class(kernel)))
      stop("An invalid kernel object is provided.")
    if (ncol(X) != kernel$get_dimension())
      stop(sprintf("Dimension of X (%d) does not match with the kernel dimension (%d).", 
                   ncol(X), kernel$get_dimension()))
  }
  # initialize kriging class
  kriging <- new(LimitKriging, X, y, kernel, interpolation)
  if (fit) {
    # set nlopt parameters 
    kriging <- Set.Kriging.NLOPT.Parameters(kriging, nlopt.parameters)
    # fit kriging
    kriging$fit(lengthscale$lower.bound, lengthscale$upper.bound)
  } else {
    # compute kriging parameters for the given lengthscale
    kriging$set_kriging_parameters()
  }
  return (kriging)
}

#' @title 
#' Universal Kriging
#' 
#' @description
#' This functions fits the universal kriging model to the data.
#' 
#' @details
#' Universal kriging permits a more general function of mean, which can be specified using \code{basis.function}. 
#' Please see Santner et al. (2003) for details.
#' 
#' For data from deterministic computer experiments, use \code{interpolation=TRUE} and will give an interpolator. 
#' For noisy data, use \code{interpolation=FALSE}, which will give an approximator of the underlying function.
#' 
#' The kernel choices are required and can be specified by 
#' (i) providing the kernel class object to \code{kernel}
#' or (ii) specifying the kernel type and other parameters in \code{kernel.parameters}. 
#' Please see examples section of \link{Fit.Kriging} for detail usages. 
#' 
#' When the lengthscale / correlation parameters are unknown, 
#' all parameters including the constant mean can be estimated via Maximum Likelihood method by setting \code{fit=TRUE}. 
#' The initial / lower bound / upper bound of the lengthscale parameters can be provided in \code{kernel.parameters}, 
#' otherwise a good initial and range would be estimated from the data. 
#' The optimization is performed via \href{https://nlopt.readthedocs.io/en/latest/}{NLopt}, 
#' a open-source library for nonlinear optimization. 
#' All gradient-free optimization methods in \href{https://nlopt.readthedocs.io/en/latest/}{NLopt} 
#' are supported and can be specified in \code{nlopt.parameters}.
#' See \code{nloptr::nloptr.print.options()} for the list of available derivative-free algorithms (prefix with NLOPT_GN or NLOPT_LN).  
#' The maximum number of optimization steps can also be defined in \code{nlopt.parameters}.
#' Please see examples section of \link{Fit.Kriging} for detail usages.
#' 
#' @param X a matrix for input (feature)
#' @param y a vector for output (target), only one-dimensional output is supported
#' @param basis.function the basis functions for specifying the prior mean
#' @param interpolation interpolation whether to interpolate, for noisy data please set \code{interpolate=FALSE}
#' @param fit whether to fit the length scale parameters from data 
#' @param kernel a kernel class object 
#' @param kernel.parameters a list of parameters required for the kernel, if no kernel class object is provided
#' @param nlopt.parameters a list of parameters required for NLopt, including choice of optimization algorithm and maximum number of evaluation
#' 
#' @return 
#' A Universal Kriging Class Object.
#' 
#' @author 
#' Chaofan Huang and V. Roshan Joseph
#' 
#' @references
#' Santner, T. J., Williams, B. J., Notz, W. I., & Williams, B. J. (2003). \emph{The design and analysis of computer experiments (Vol. 1)}. New York: Springer.
#' 
#' @seealso 
#' \link{Fit.Kriging}, \link{Predict.Kriging}, \link{Get.Kriging.Parameters}.
#' 
#' @export
#' 
#' @examples
#' # one dimensional example 
#' f <- function(x) {
#'   x <- 0.5 + 2*x
#'   y <- sin(10*pi*x)/(2*x) + (x-1)^4
#'   return (y)
#' }
#' 
#' set.seed(1234)
#' # train set
#' n <- 30
#' p <- 1
#' X <- matrix(runif(n),ncol=p)
#' y <- apply(X, 1, f)
#  # test set
#' newX <- matrix(seq(0,1,length=1001), ncol=p)
#' 
#' basis.function <- function(x) {c(1,x[1],x[1]^2)}
#' 
#' # approach 1
#' kriging <- Universal.Kriging(X, y, basis.function=basis.function,
#'                              interpolation=TRUE, fit=TRUE, 
#'                              kernel.parameters=list(type="Gaussian"))
#' pred <- Predict.Kriging(kriging, newX)
#' plot(newX, f(newX), "l")
#' points(X, y, pch=16, col="blue")
#' lines(newX, pred$mean, col="red", lty=2)
#' lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#' lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#' Get.Kriging.Parameters(kriging)
#' 
#' # approach 2
#' kriging <- Fit.Kriging(X, y, interpolation=TRUE, fit=TRUE, model="UK",
#'                        model.parameters=list(basis.function=basis.function),
#'                        kernel.parameters=list(type="Gaussian"))
#' pred <- Predict.Kriging(kriging, newX)
#' plot(newX, f(newX), "l")
#' points(X, y, pch=16, col="blue")
#' lines(newX, pred$mean, col="red", lty=2)
#' lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#' lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#' Get.Kriging.Parameters(kriging)
Universal.Kriging <- function(X, y, basis.function, interpolation=TRUE, fit=TRUE,
                              kernel=NULL, kernel.parameters=list(), 
                              nlopt.parameters=list()) {
  # check inputs
  if (!is.matrix(X)) X <- as.matrix(X)
  if (is.matrix(y)) y <- c(y)
  # get lengthscale parameters if not provided
  lengthscale <- list(
    value = kernel.parameters[["lengthscale"]],
    lower.bound = kernel.parameters[["lengthscale.lower.bound"]],
    upper.bound = kernel.parameters[["lengthscale.upper.bound"]]
  )
  if (fit) lengthscale <- Get.Lengthscale(X, lengthscale)
  # get kernel class
  if (is.null(kernel)) {
    if (!("type" %in% names(kernel.parameters))) 
      stop("either kernel object or type of kernel must be provided in kernel.parameters.")
    if (!("value" %in% names(lengthscale)))
      stop("either kernel object or lengthscale of kernel must be provided in kernel.parameters when fit = FALSE.")
    kernel <- Get.Kernel(lengthscale$value, kernel.parameters[["type"]], kernel.parameters)
  } else {
    if (!grepl("^Rcpp.*Kernel$", class(kernel)))
      stop("An invalid kernel object is provided.")
    if (ncol(X) != kernel$get_dimension())
      stop(sprintf("Dimension of X (%d) does not match with the kernel dimension (%d).", 
                   ncol(X), kernel$get_dimension()))
  }
  # initialize kriging class
  if (!is.function(basis.function))
    stop("basis.function must a function for Universal Kriging (UK).")
  no.basis.function <- length(basis.function(X[1,]))
  kriging <- new(UniversalKriging, X, y, kernel, interpolation, no.basis.function, basis.function)
  if (fit) {
    # set nlopt parameters 
    kriging <- Set.Kriging.NLOPT.Parameters(kriging, nlopt.parameters)
    # fit kriging
    kriging$fit(lengthscale$lower.bound, lengthscale$upper.bound)
  } else {
    # compute kriging parameters for the given lengthscale
    kriging$set_kriging_parameters()
  }
  return (kriging)
}

#' @title 
#' Rational Kriging
#' 
#' @description
#' This functions fits the rational kriging model to the data.
#' 
#' @details
#' Rational kriging gives a rational predictor in terms of the inputs, which gives a more stable estimate of the mean. 
#' Please see Joseph (2024) for details.
#' Only interpolation is available. Noisy output is not supported for rational kriging.  
#' 
#' The kernel choices are required and can be specified by 
#' (i) providing the kernel class object to \code{kernel}
#' or (ii) specifying the kernel type and other parameters in \code{kernel.parameters}. 
#' Please see examples section of \link{Fit.Kriging} for detail usages. 
#' 
#' When the lengthscale / correlation parameters are unknown, 
#' all parameters including the constant mean can be estimated via Maximum Likelihood method by setting \code{fit=TRUE}. 
#' The initial / lower bound / upper bound of the lengthscale parameters can be provided in \code{kernel.parameters}, 
#' otherwise a good initial and range would be estimated from the data. 
#' The optimization is performed via \href{https://nlopt.readthedocs.io/en/latest/}{NLopt}, 
#' a open-source library for nonlinear optimization. 
#' All gradient-free optimization methods in \href{https://nlopt.readthedocs.io/en/latest/}{NLopt} 
#' are supported and can be specified in \code{nlopt.parameters}.
#' See \code{nloptr::nloptr.print.options()} for the list of available derivative-free algorithms (prefix with NLOPT_GN or NLOPT_LN).  
#' The maximum number of optimization steps can also be defined in \code{nlopt.parameters}.
#' Please see examples section of \link{Fit.Kriging} for detail usages.
#' 
#' @param X a matrix for input (feature)
#' @param y a vector for output (target), only one-dimensional output is supported
#' @param fit whether to fit the length scale parameters from data 
#' @param kernel a kernel class object 
#' @param kernel.parameters a list of parameters required for the kernel, if no kernel class object is provided
#' @param nlopt.parameters a list of parameters required for NLopt, including choice of optimization algorithm and maximum number of evaluation
#' 
#' @return 
#' A Rational Kriging Class Object.
#' 
#' @author 
#' Chaofan Huang and V. Roshan Joseph
#' 
#' @references
#' Joseph, V. R. (2024). Rational Kriging. \emph{Journal of the American Statistical Association}.
#' 
#' @seealso 
#' \link{Fit.Kriging}, \link{Predict.Kriging}, \link{Get.Kriging.Parameters}.
#' 
#' @export
#' 
#' @examples
#' # one dimensional example 
#' f <- function(x) {
#'   x <- 0.5 + 2*x
#'   y <- sin(10*pi*x)/(2*x) + (x-1)^4
#'   return (y)
#' }
#' 
#' set.seed(1234)
#' # train set
#' n <- 30
#' p <- 1
#' X <- matrix(runif(n),ncol=p)
#' y <- apply(X, 1, f)
#  # test set
#' newX <- matrix(seq(0,1,length=1001), ncol=p)
#' 
#' # approach 1
#' kriging <- Rational.Kriging(X, y, fit=TRUE, kernel.parameters=list(type="RQ",alpha=1))
#' pred <- Predict.Kriging(kriging, newX)
#' plot(newX, f(newX), "l")
#' points(X, y, pch=16, col="blue")
#' lines(newX, pred$mean, col="red", lty=2)
#' lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#' lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#' Get.Kriging.Parameters(kriging)
#' 
#' # approach 2
#' kriging <- Fit.Kriging(X, y, interpolation=TRUE, fit=TRUE, model="RK",
#'                        kernel.parameters=list(type="RQ",alpha=1))
#' pred <- Predict.Kriging(kriging, newX)
#' plot(newX, f(newX), "l")
#' points(X, y, pch=16, col="blue")
#' lines(newX, pred$mean, col="red", lty=2)
#' lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#' lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#' Get.Kriging.Parameters(kriging)
#' 
Rational.Kriging <- function(X, y, fit=TRUE,
                             kernel=NULL, kernel.parameters=list(), 
                             nlopt.parameters=list()) {
  # only interpolation is supported for the first release
  interpolation <- TRUE
  # check inputs
  if (!is.matrix(X)) X <- as.matrix(X)
  if (is.matrix(y)) y <- c(y)
  # get lengthscale parameters if not provided
  lengthscale <- list(
    value = kernel.parameters[["lengthscale"]],
    lower.bound = kernel.parameters[["lengthscale.lower.bound"]],
    upper.bound = kernel.parameters[["lengthscale.upper.bound"]]
  )
  if (fit) lengthscale <- Get.Lengthscale(X, lengthscale)
  # get kernel class
  if (is.null(kernel)) {
    if (!("type" %in% names(kernel.parameters))) 
      stop("either kernel object or type of kernel must be provided in kernel.parameters.")
    if (!("value" %in% names(lengthscale)))
      stop("either kernel object or lengthscale of kernel must be provided in kernel.parameters when fit = FALSE.")
    kernel <- Get.Kernel(lengthscale$value, kernel.parameters[["type"]], kernel.parameters)
  } else {
    if (!grepl("^Rcpp.*Kernel$", class(kernel)))
      stop("An invalid kernel object is provided.")
    if (ncol(X) != kernel$get_dimension())
      stop(sprintf("Dimension of X (%d) does not match with the kernel dimension (%d).", 
                   ncol(X), kernel$get_dimension()))
  }
  # initialize kriging class
  kriging <- new(RationalKriging, X, y, kernel, interpolation)
  if (fit) {
    # set nlopt parameters 
    kriging <- Set.Kriging.NLOPT.Parameters(kriging, nlopt.parameters)
    # fit kriging
    kriging$fit(lengthscale$lower.bound, lengthscale$upper.bound)
  } else {
    # compute kriging parameters for the given lengthscale
    kriging$set_kriging_parameters()
  }
  return (kriging)
}

#' @title 
#' Generalized Rational Kriging
#' 
#' @description
#' This functions fits the generalized rational kriging model to the data.
#' 
#' @details
#' Ordinary kriging and rational kriging can be obtained as special cases of generalized rational kriging. 
#' Please see Joseph (2024) for details. 
#' The \href{https://spectralib.org/}{Spectra} library is used for fast computation of the first eigenvalues/vectors. 
#' Only interpolation is available. Noisy output is not supported for generalized rational kriging.  
#' 
#' The kernel choices are required and can be specified by 
#' (i) providing the kernel class object to \code{kernel}
#' or (ii) specifying the kernel type and other parameters in \code{kernel.parameters}. 
#' Please see examples section of \link{Fit.Kriging} for detail usages. 
#' 
#' When the lengthscale / correlation parameters are unknown, 
#' all parameters including the constant mean can be estimated via Maximum Likelihood method by setting \code{fit=TRUE}. 
#' The initial / lower bound / upper bound of the lengthscale parameters can be provided in \code{kernel.parameters}, 
#' otherwise a good initial and range would be estimated from the data. 
#' The optimization is performed via \href{https://nlopt.readthedocs.io/en/latest/}{NLopt} library, 
#' a open-source library for nonlinear optimization. 
#' All gradient-free optimization methods in \href{https://nlopt.readthedocs.io/en/latest/}{NLopt} 
#' are supported and can be specified in \code{nlopt.parameters}.
#' See \code{nloptr::nloptr.print.options()} for the list of available derivative-free algorithms (prefix with NLOPT_GN or NLOPT_LN).  
#' The maximum number of optimization steps can also be defined in \code{nlopt.parameters}.
#' Please see examples section of \link{Fit.Kriging} for detail usages.
#' 
#' @param X a matrix for input (feature)
#' @param y a vector for output (target), only one-dimensional output is supported
#' @param fit whether to fit the length scale parameters from data 
#' @param kernel a kernel class object 
#' @param kernel.parameters a list of parameters required for the kernel, if no kernel class object is provided
#' @param nlopt.parameters a list of parameters required for NLopt, including choice of optimization algorithm and maximum number of evaluation
#' 
#' @return 
#' A Generalized Rational Kriging Class Object.
#' 
#' @author 
#' Chaofan Huang and V. Roshan Joseph
#' 
#' @references
#' Joseph, V. R. (2024). Rational Kriging. \emph{Journal of the American Statistical Association}.
#' 
#' Qiu, Y., Guennebaud, G., & Niesen, J. (2015). Spectra: C++ library for large scale eigenvalue problems.
#' 
#' @seealso 
#' \link{Fit.Kriging}, \link{Predict.Kriging}, \link{Get.Kriging.Parameters}.
#' 
#' @export
#' 
#' @examples
#' # one dimensional example 
#' f <- function(x) {
#'   x <- 0.5 + 2*x
#'   y <- sin(10*pi*x)/(2*x) + (x-1)^4
#'   return (y)
#' }
#' 
#' set.seed(1234)
#' # train set
#' n <- 30
#' p <- 1
#' X <- matrix(runif(n),ncol=p)
#' y <- apply(X, 1, f)
#  # test set
#' newX <- matrix(seq(0,1,length=1001), ncol=p)
#' 
#' # approach 1
#' kriging <- Generalized.Rational.Kriging(X, y, fit=TRUE, 
#'                                         kernel.parameters=list(type="RQ",alpha=1))
#' pred <- Predict.Kriging(kriging, newX)
#' plot(newX, f(newX), "l")
#' points(X, y, pch=16, col="blue")
#' lines(newX, pred$mean, col="red", lty=2)
#' lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#' lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#' Get.Kriging.Parameters(kriging)
#' 
#' # approach 2
#' kriging <- Fit.Kriging(X, y, interpolation=TRUE, fit=TRUE, model="GRK",
#'                        kernel.parameters=list(type="RQ",alpha=1))
#' pred <- Predict.Kriging(kriging, newX)
#' plot(newX, f(newX), "l")
#' points(X, y, pch=16, col="blue")
#' lines(newX, pred$mean, col="red", lty=2)
#' lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#' lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#' Get.Kriging.Parameters(kriging)
#' 
Generalized.Rational.Kriging <- function(X, y, fit=TRUE,  
                                         kernel=NULL, kernel.parameters=list(), 
                                         nlopt.parameters=list()) {
  # only interpolation is supported for the first release
  interpolation <- TRUE
  # check inputs
  if (!is.matrix(X)) X <- as.matrix(X)
  if (is.matrix(y)) y <- c(y)
  # get lengthscale parameters if not provided
  lengthscale <- list(
    value = kernel.parameters[["lengthscale"]],
    lower.bound = kernel.parameters[["lengthscale.lower.bound"]],
    upper.bound = kernel.parameters[["lengthscale.upper.bound"]]
  )
  if (fit) lengthscale <- Get.Lengthscale(X, lengthscale)
  # get kernel class
  if (is.null(kernel)) {
    if (!("type" %in% names(kernel.parameters))) 
      stop("either kernel object or type of kernel must be provided in kernel.parameters.")
    if (!("value" %in% names(lengthscale)))
      stop("either kernel object or lengthscale of kernel must be provided in kernel.parameters when fit = FALSE.")
    kernel <- Get.Kernel(lengthscale$value, kernel.parameters[["type"]], kernel.parameters)
  } else {
    if (!grepl("^Rcpp.*Kernel$", class(kernel)))
      stop("An invalid kernel object is provided.")
    if (ncol(X) != kernel$get_dimension())
      stop(sprintf("Dimension of X (%d) does not match with the kernel dimension (%d).", 
                   ncol(X), kernel$get_dimension()))
  }
  # initialize kriging class
  kriging <- new(GeneralizedRationalKriging, X, y, kernel, interpolation)
  if (fit) {
    # set nlopt parameters 
    kriging <- Set.Kriging.NLOPT.Parameters(kriging, nlopt.parameters)
    # fit kriging
    kriging$fit(lengthscale$lower.bound, lengthscale$upper.bound)
  } else {
    # compute kriging parameters for the given lengthscale
    kriging$set_kriging_parameters()
  }
  return (kriging)
}

#' @title 
#' Kriging Prediction
#' 
#' @description
#' This function gives prediction and uncertainty quantification of the kriging model on a new input.
#' 
#' @param kriging a kriging class object
#' @param X a matrix for the new input (features) to perform predictions
#' 
#' @return 
#' \item{mean}{kriging mean computed at the new input}
#' \item{sd}{kriging standard computed at the new input}
#' 
#' @author 
#' Chaofan Huang and V. Roshan Joseph
#' 
#' @references
#' Joseph, V. R. (2006). \emph{Limit kriging}. Technometrics, 48(4), 458-466.
#' 
#' Joseph, V. R. (2024). Rational Kriging. \emph{Journal of the American Statistical Association}.
#' 
#' Rasmussen, C. E. & Williams, C. K. (2006). \emph{Gaussian Processes for Machine Learning}. The MIT Press.
#' 
#' Santner, T. J., Williams, B. J., Notz, W. I., & Williams, B. J. (2003). \emph{The design and analysis of computer experiments (Vol. 1)}. New York: Springer.
#' 
#' @seealso 
#' \link{Fit.Kriging}.
#' 
#' @export
#' 
#' @examples 
#' # one dimensional example 
#' f <- function(x) {
#'   x <- 0.5 + 2*x
#'   y <- sin(10*pi*x)/(2*x) + (x-1)^4
#'   return (y)
#' }
#' 
#' set.seed(1234)
#' # train set
#' n <- 30
#' p <- 1
#' X <- matrix(runif(n),ncol=p)
#' y <- apply(X, 1, f)
#  # test set
#' newX <- matrix(seq(0,1,length=1001), ncol=p)
#' 
#' kriging <- Fit.Kriging(X, y, interpolation=TRUE, fit=TRUE, model="OK",
#'                        kernel.parameters=list(type="Gaussian"))
#' pred <- Predict.Kriging(kriging, newX)
#' plot(newX, f(newX), "l")
#' points(X, y, pch=16, col="blue")
#' lines(newX, pred$mean, col="red", lty=2)
#' lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#' lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#' 
Predict.Kriging <- function(kriging, X) {
  if (!grepl("^Rcpp.*Kriging$", class(kriging)))
    stop("An invalid kriging object is provided.")
  if (!is.matrix(X)) X <- as.matrix(X)
  if (ncol(X) != kriging$get_dimension())
    stop(sprintf("Dimension of X (%d) does not match with the kriging dimension (%d).", 
                 ncol(X), kriging$get_dimension()))
  return (kriging$predict(X))
}

#' @title 
#' Get Kriging Parameters
#' 
#' @description
#' This function can be used for extracting the estimates of the kriging parameters.
#' 
#' @param kriging a kriging class object
#' 
#' @return 
#' \item{nllh}{negative log-likelihood of the kriging model}
#' \item{mu}{mean of the kriging model}
#' \item{nu2}{variance of the kriging model}
#' \item{sigma2}{variance of the random noise when \code{interpolation=FALSE}}
#' \item{beta}{coefficients of the basis functions for universal kriging}
#' \item{c}{c for the rational / generalized rational kriging, see Joseph (2024)}
#' \item{c0}{c0 for the generalized rational kriging, see Joseph (2024)}
#' 
#' @author 
#' Chaofan Huang and V. Roshan Joseph
#' 
#' @seealso 
#' \link{Fit.Kriging}.
#' 
#' @references
#' Joseph, V. R. (2006). \emph{Limit kriging}. Technometrics, 48(4), 458-466.
#' 
#' Joseph, V. R. (2024). Rational Kriging. \emph{Journal of the American Statistical Association}.
#' 
#' Rasmussen, C. E. & Williams, C. K. (2006). \emph{Gaussian Processes for Machine Learning}. The MIT Press.
#' 
#' Santner, T. J., Williams, B. J., Notz, W. I., & Williams, B. J. (2003). \emph{The design and analysis of computer experiments (Vol. 1)}. New York: Springer.
#' 
#' @export
#' 
#' @examples 
#' # one dimensional example 
#' f <- function(x) {
#'   x <- 0.5 + 2*x
#'   y <- sin(10*pi*x)/(2*x) + (x-1)^4
#'   return (y)
#' }
#' 
#' set.seed(1234)
#' # train set
#' n <- 30
#' p <- 1
#' X <- matrix(runif(n),ncol=p)
#' y <- apply(X, 1, f)
#  # test set
#' newX <- matrix(seq(0,1,length=1001), ncol=p)
#' 
#' kriging <- Fit.Kriging(X, y, interpolation=TRUE, fit=TRUE, model="OK",
#'                        kernel.parameters=list(type="Gaussian"))
#' Get.Kriging.Parameters(kriging)
#' 
Get.Kriging.Parameters <- function(kriging) {
  if (!grepl("^Rcpp.*Kriging$", class(kriging)))
    stop("An invalid kriging object is provided.")
  kriging.parameters <- list(
    "nllh" = kriging$get_nllh(),
    "lengthscale" = kriging$get_lengthscale(),
    "mu" = kriging$get_mu(),
    "nu2" = kriging$get_nu2(),
    "sigma2" = kriging$get_sigma2()
  )
  if (inherits(kriging, "Rcpp_UniversalKriging")) {
    kriging.parameters[["beta"]] <- kriging$get_beta()
    kriging.parameters[["mu"]] <- NULL
  }
  if (inherits(kriging, "Rcpp_RationalKriging"))
    kriging.parameters[["c"]] <- kriging$get_c()
  if (inherits(kriging, "Rcpp_GeneralizedRationalKriging")) {
    kriging.parameters[["c0"]] <- kriging$get_c0()
    kriging.parameters[["c"]] <- kriging$get_c()
  }
  return (kriging.parameters)
}
