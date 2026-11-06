#' Sequential Class
#'
#' @useDynLib Sequential, .registration = TRUE
#' @import Rcpp, RcppEigen, BH, nloptr
#' 
#' @noRd
Rcpp::loadModule("Sequential", TRUE)

#' @noRd
Set.Sequential.NLOPT.Parameters <- function(sequential, nlopt.parameters=list()) {
  if (!grepl("^Rcpp.*Sequential$", class(sequential)))
    stop("An invalid sequential object is provided.")
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
    sequential$set_nlopt_algorithm(nlopt.algorithm)
    if (nlopt.algorithm %in% c("NLOPT_GN_MLSL","NLOPT_GN_MLSL_LDS")) {
      if ("local.algorithm" %in% names(nlopt.parameters)) {
        nlopt.local.algorithm <- nlopt.parameters[["local.algorithm"]]
        if (!(nlopt.local.algorithm %in% local.optimization.algorithms))
          stop(sprintf("%s optimization algorithm is not supported. Only local gradient-free optimization is supported: %s", 
                       nlopt.local.algorithm, paste(local.optimization.algorithms, collapse=", ")))
        sequential$set_nlopt_local_algorithm(nlopt.local.algorithm)
      } else {
        sequential$set_nlopt_local_algorithm("NLOPT_LN_SBPLX")
      }
    }
  } 
  if ("maxeval" %in% names(nlopt.parameters)) {
    nlopt.maxeval <- as.integer(nlopt.parameters[["maxeval"]])
    if (nlopt.maxeval <= 0) 
      stop(sprintf("nlopt optimization maximum evaluation must be a positive integer."))
    sequential$set_nlopt_maxeval(nlopt.maxeval)
  }
  if ("population" %in% names(nlopt.parameters)) {
    nlopt.population <- as.integer(nlopt.parameters[["population"]])
    if (nlopt.population <= 0) 
      stop(sprintf("nlopt optimization population size must be a positive integer."))
    sequential$set_nlopt_population(nlopt.population)
  }
  return (sequential)
}

#' @noRd
Get.Sequential.Design.Setup <- function(design.parameters, kriging) {
  if (!is.list(design.parameters))
    stop("design.parameters must be a list.")
  if (!all(c("design.lower.bound", "design.upper.bound") %in% names(design.parameters)))
    stop("design.parameters must contain design.lower.bound and design.upper.bound.")
  p <- kriging$get_dimension()
  design.lower.bound <- c(design.parameters[["design.lower.bound"]])
  design.upper.bound <- c(design.parameters[["design.upper.bound"]])
  if (!is.numeric(design.lower.bound) || !is.numeric(design.upper.bound) ||
      length(design.lower.bound) != p || length(design.upper.bound) != p ||
      any(!is.finite(design.lower.bound)) || any(!is.finite(design.upper.bound)))
    stop(sprintf("Design bounds must be finite numeric vectors matching the kriging dimension (%d).", p))
  if (any(design.lower.bound > design.upper.bound))
    stop("design.lower.bound can not be greater than design.upper.bound.")
  scale.info <- Get.Scale.Info(kriging, p, "kriging")
  if (!is.null(scale.info)) {
    design.lower.bound <- c(Scale.Kriging.Input(matrix(design.lower.bound, nrow=1), scale.info$min, scale.info$max))
    design.upper.bound <- c(Scale.Kriging.Input(matrix(design.upper.bound, nrow=1), scale.info$min, scale.info$max))
  }
  return (list(
    lower.bound=design.lower.bound,
    upper.bound=design.upper.bound,
    scale.info=scale.info
  ))
}

#' @noRd
Set.Sequential.Scale.Info <- function(sequential, scale.info=NULL) {
  attr(sequential, "scale") <- !is.null(scale.info)
  if (!is.null(scale.info)) {
    attr(sequential, "scale_min") <- scale.info$min
    attr(sequential, "scale_max") <- scale.info$max
  }
  return (sequential)
}

#' @noRd
Scale.Sequential.Input <- function(sequential, X) {
  scale.info <- Get.Scale.Info(sequential, ncol(X), "sequential")
  if (is.null(scale.info))
    return (X)
  if (!is.numeric(X) || any(!is.finite(X)))
    stop("X must be a finite numeric matrix for a scaled sequential object.")
  return (Scale.Kriging.Input(X, scale.info$min, scale.info$max))
}

#' @noRd
Unscale.Sequential.Input <- function(sequential, X) {
  scale.info <- Get.Scale.Info(sequential, ncol(X), "sequential")
  if (is.null(scale.info))
    return (X)
  return (Unscale.Kriging.Input(X, scale.info$min, scale.info$max))
}

#' @noRd
Wrap.Scaled.Sequential.Criterion <- function(criterion.function, scale.info=NULL) {
  if (is.null(scale.info))
    return (criterion.function)
  force(criterion.function)
  force(scale.info)
  function(mean, sd, X, y, criterion.parameters) {
    X <- Unscale.Kriging.Input(X, scale.info$min, scale.info$max)
    criterion.function(mean, sd, X, y, criterion.parameters)
  }
}

#' @title
#' Get Sequential Design
#'
#' @description
#' This function provides a common interface to create various sequential design objects.
#' ALM (\link{ALM.Sequential}), EI (\link{EI.Sequential}), PI (\link{PI.Sequential}),
#' UCB (\link{UCB.Sequential}), and UDF (\link{UDF.Sequential}) are supported.
#'
#' @details
#' Sequential design methods select the next design point iteratively based on an acquisition criterion.
#' For emulation tasks (space-filling), use \code{method="ALM"} which maximizes predictive uncertainty.
#' For optimization tasks (finding the minimum), use \code{method="EI"}, \code{"PI"}, or \code{"UCB"}.
#' A custom criterion can be specified via \code{method="UDF"}.
#'
#' The design space is specified via \code{design.parameters}, which must contain
#' \code{design.lower.bound} and \code{design.upper.bound} vectors.
#' During optimization, candidate points can be provided via \code{design.candidates}
#' in \code{\link{Get.Next.Design}}.
#' If the kriging object was fitted with \code{scale=TRUE}, design bounds,
#' candidates, new design points, predictions, returned designs, and returned data
#' are supplied and returned on the original input scale; scaling is handled
#' internally for the sequential criterion.
#' For \code{method="UDF"}, the user-defined criterion function receives the
#' current training \code{X} on the original input scale.
#'
#' The kriging model can be provided either as:
#' \describe{
#'   \item{An existing kriging object}{e.g., from \link{Ordinary.Kriging} or \link{Fit.Kriging}}
#'   \item{Kriging parameters}{as a list to pass to \link{Fit.Kriging}, must contain at least X and y}
#' }
#'
#' @param design.parameters a list specifying the design space, containing:
#' \describe{
#'   \item{\code{design.lower.bound}}{a vector for the lower bound of the design space}
#'   \item{\code{design.upper.bound}}{a vector for the upper bound of the design space}
#' }
#' @param kriging (optional if kriging.parameters provided) a fitted kriging class object 
#' @param kriging.parameters (optional if kriging provided) a list of parameters to create kriging via \link{Fit.Kriging},
#'        must contain X, y and optionally other \link{Fit.Kriging} parameters
#' @param method choice of sequential design method: ALM, EI, PI, UCB, UDF
#' @param method.parameters a list of parameters required for the specific method, e.g. \code{lambda} for UCB,
#'        \code{criterion.function} and \code{criterion.parameters} for UDF
#' @param nlopt.parameters a list of parameters for NLopt optimization (used when \code{design.candidates} is not provided)
#'
#' @return
#' A Sequential Class Object.
#'
#' @author
#' Chaofan Huang, Difan Song, and V. Roshan Joseph
#'
#' @seealso
#' \link{ALM.Sequential}, \link{EI.Sequential}, \link{PI.Sequential}, \link{UCB.Sequential},
#' \link{UDF.Sequential}, \link{Predict.Sequential},
#' \link{Get.Next.Design}, \link{Add.Next.Design}.
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
#' n <- 10
#' p <- 1
#' X <- matrix(runif(n), ncol=p)
#' y <- apply(X, 1, f)
#' newX <- matrix(seq(0, 1, length=1001), ncol=p)
#'
#' # create sequential
#' kriging <- Ordinary.Kriging(X, y, interpolation=TRUE, fit=TRUE,
#'                             kernel.parameters=list(type="Gaussian"))
#' design.parameters <- list(design.lower.bound=c(0), design.upper.bound=c(1))
#' sequential <- Get.Sequential(design.parameters, kriging, method="EI")
#'
#' # sequential design loop with visualization
#' for (i in 1:5) {
#'   # visualize current kriging
#'   pred <- Predict.Sequential(sequential, newX)
#'   plot(newX, f(newX), "l", main=sprintf("EI Iteration %d", i))
#'   points(Get.Data(sequential)$X, Get.Data(sequential)$y, pch=16, col="blue")
#'   lines(newX, pred$mean, col="red", lty=2)
#'   lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#'   lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#'   # get next design and evaluate
#'   Xn <- Get.Next.Design(sequential)
#'   new.design <- list(X=Xn, y=f(Xn))
#'   points(new.design$X, new.design$y, pch=16, col="green")
#'   # add to kriging and update
#'   sequential <- Add.Next.Design(sequential, new.design)
#' }
#'
#' # final visualization
#' pred <- Predict.Sequential(sequential, newX)
#' plot(newX, f(newX), "l", main="EI Final")
#' points(Get.Data(sequential)$X, Get.Data(sequential)$y, pch=16, col="blue")
#' lines(newX, pred$mean, col="red", lty=2)
#' lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#' lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#' best.design <- Get.Best.Design(sequential)
#' points(best.design$design, best.design$objective, pch=16, col="green")
#'
#' # alternatively can create sequential by passing kriging parameters
#' design.parameters <- list(design.lower.bound=c(0), design.upper.bound=c(1))
#' kriging.parameters <- list(X=X, y=y, interpolation=TRUE, fit=TRUE,
#'                            kernel.parameters=list(type="Gaussian"))
#' sequential <- Get.Sequential(design.parameters, kriging.parameters=kriging.parameters, method="EI")
#' 
Get.Sequential <- function(design.parameters, kriging=NULL, kriging.parameters=NULL,
                           method="ALM", method.parameters=list(),
                           nlopt.parameters=list()) {
  # check that at least one of kriging or kriging.parameters is provided
  if (is.null(kriging) & is.null(kriging.parameters))
    stop("Either kriging object or kriging.parameters must be provided.")
  
  # create kriging if kriging is not provided
  if (is.null(kriging)) {
    if (!is.list(kriging.parameters))
      stop("kriging.parameters must be a list.")
    if (!all(c("X", "y") %in% names(kriging.parameters)))
      stop("kriging.parameters must contain X and y.")
    # Call Fit.Kriging with the provided parameters
    kriging <- do.call(Fit.Kriging, kriging.parameters)
  }
  
  # validate kriging object
  if (!grepl("^Rcpp.*Kriging$", class(kriging)))
    stop("An invalid kriging object is provided.")
  
  design.setup <- Get.Sequential.Design.Setup(design.parameters, kriging)
  
  # create sequential based on method
  if (!is.character(method))
    stop("method must be a string.")
  if (method == "ALM") {
    sequential <- new(ALMSequential, kriging, design.setup$lower.bound, design.setup$upper.bound)
  } else if (method == "EI") {
    sequential <- new(EISequential, kriging, design.setup$lower.bound, design.setup$upper.bound)
  } else if (method == "PI") {
    sequential <- new(PISequential, kriging, design.setup$lower.bound, design.setup$upper.bound)
  } else if (method == "UCB") {
    if (!("lambda" %in% names(method.parameters)))
      stop("lambda must be provided in method.parameters for UCB method.")
    lambda <- method.parameters[["lambda"]]
    if (!(is.numeric(lambda) & lambda > 0))
      stop("lambda must be a positive number.")
    sequential <- new(UCBSequential, kriging, design.setup$lower.bound, design.setup$upper.bound, lambda)
  } else if (method == "UDF") {
    if (!("criterion.function" %in% names(method.parameters)))
      stop("criterion.function must be provided in method.parameters for UDF method.")
    if (!is.function(method.parameters[["criterion.function"]]))
      stop("criterion.function must be an R function.")
    criterion.parameters <- if ("criterion.parameters" %in% names(method.parameters)) method.parameters[["criterion.parameters"]] else list()
    if (!is.list(criterion.parameters))
      stop("criterion.parameters must be a list.")
    criterion.function <- Wrap.Scaled.Sequential.Criterion(method.parameters[["criterion.function"]], design.setup$scale.info)
    sequential <- new(UDFSequential, kriging, design.setup$lower.bound, design.setup$upper.bound, criterion.function, criterion.parameters)
  } else {
    stop(sprintf("No sequential method '%s'. Supported methods: ALM, EI, PI, UCB, UDF.", method))
  }
  # set nlopt parameters
  sequential <- Set.Sequential.Scale.Info(sequential, design.setup$scale.info)
  sequential <- Set.Sequential.NLOPT.Parameters(sequential, nlopt.parameters)
  return (sequential)
}

#' @title
#' ALM Sequential Design
#' 
#' @description
#' This function creates an Active Learning MacKay (ALM) sequential design object 
#' for emulation (space-filling based on predictive uncertainty).
#' 
#' @details
#' ALM selects the next design point as the one that maximizes the predictive standard deviation:
#' \deqn{x_{next} = \arg\max_{x} \sigma(x),}{x_next = argmax_x sigma(x),}
#' where \eqn{\sigma(x)}{sigma(x)} is the kriging predictive standard deviation at \eqn{x}{x}.
#' This criterion is suitable for emulation tasks where the goal is to reduce overall prediction uncertainty.
#'
#' The design space is specified via \code{design.parameters}, which must contain
#' \code{design.lower.bound} and \code{design.upper.bound} vectors.
#' During optimization, candidate points can be provided via \code{design.candidates}
#' in \code{\link{Get.Next.Design}}.
#' If the kriging object was fitted with \code{scale=TRUE}, design bounds,
#' candidates, new design points, predictions, returned designs, and returned data
#' are supplied and returned on the original input scale; scaling is handled
#' internally for the sequential criterion.
#'
#' @param design.parameters a list specifying the design space, containing:
#' \describe{
#'   \item{\code{design.lower.bound}}{a vector for the lower bound of the design space}
#'   \item{\code{design.upper.bound}}{a vector for the upper bound of the design space}
#' }
#' @param kriging (optional if kriging.parameters provided) a fitted kriging class object (e.g., from \link{Ordinary.Kriging} or \link{Fit.Kriging})
#' @param kriging.parameters (optional if kriging provided) a list of parameters to create kriging via \link{Fit.Kriging},
#'        must contain X, y and optionally other \link{Fit.Kriging} parameters
#' @param nlopt.parameters a list of parameters for NLopt optimization (used when \code{design.candidates} is not provided)
#' 
#' @return 
#' An ALM Sequential Class Object.
#' 
#' @author 
#' Chaofan Huang, Difan Song, and V. Roshan Joseph
#' 
#' @references
#' MacKay, D. J. C. (1992). Information-based objective functions for active data selection. 
#' \emph{Neural computation}, 4(4), 590-604.
#' 
#' @seealso 
#' \link{Predict.Sequential}, \link{Get.Next.Design}, \link{Add.Next.Design}.
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
#' n <- 10
#' p <- 1
#' X <- matrix(runif(n), ncol=p)
#' y <- apply(X, 1, f)
#' newX <- matrix(seq(0, 1, length=1001), ncol=p)
#'
#' # create sequential
#' kriging <- Ordinary.Kriging(X, y, interpolation=TRUE, fit=TRUE,
#'                             kernel.parameters=list(type="Gaussian"))
#' design.parameters <- list(design.lower.bound=c(0), design.upper.bound=c(1))
#' sequential <- ALM.Sequential(design.parameters, kriging)
#'
#' # sequential design loop with visualization
#' n.iterations <- 5
#' for (i in 1:n.iterations) {
#'   # visualize current kriging
#'   pred <- Predict.Sequential(sequential, newX)
#'   plot(newX, f(newX), "l", main=sprintf("ALM Iteration %d", i))
#'   points(Get.Data(sequential)$X, Get.Data(sequential)$y, pch=16, col="blue")
#'   lines(newX, pred$mean, col="red", lty=2)
#'   lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#'   lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#'   # get next design and evaluate
#'   Xn <- Get.Next.Design(sequential)
#'   new.design <- list(X=Xn, y=f(Xn))
#'   points(new.design$X, new.design$y, pch=16, col="green")
#'   # add to kriging and update
#'   sequential <- Add.Next.Design(sequential, new.design)
#' }
#'
#' # final visualization
#' pred <- Predict.Sequential(sequential, newX)
#' plot(newX, f(newX), "l", main="ALM Final")
#' points(Get.Data(sequential)$X, Get.Data(sequential)$y, pch=16, col="blue")
#' lines(newX, pred$mean, col="red", lty=2)
#' lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#' lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#'
ALM.Sequential <- function(design.parameters, kriging=NULL, kriging.parameters=NULL, nlopt.parameters=list()) {
  # check that at least one of kriging or kriging.parameters is provided
  if (is.null(kriging) & is.null(kriging.parameters))
    stop("Either kriging object or kriging.parameters must be provided.")
  
  # create kriging if kriging is not provided
  if (is.null(kriging)) {
    if (!is.list(kriging.parameters))
      stop("kriging.parameters must be a list.")
    if (!all(c("X", "y") %in% names(kriging.parameters)))
      stop("kriging.parameters must contain X and y.")
    kriging <- do.call(Fit.Kriging, kriging.parameters)
  }
  
  if (!grepl("^Rcpp.*Kriging$", class(kriging)))
    stop("An invalid kriging object is provided.")
  design.setup <- Get.Sequential.Design.Setup(design.parameters, kriging)
  sequential <- new(ALMSequential, kriging, design.setup$lower.bound, design.setup$upper.bound)
  sequential <- Set.Sequential.Scale.Info(sequential, design.setup$scale.info)
  sequential <- Set.Sequential.NLOPT.Parameters(sequential, nlopt.parameters)
  return (sequential)
}

#' @title 
#' EI Sequential Design
#' 
#' @description
#' This function creates an Expected Improvement (EI) sequential design object 
#' for optimization (finding the minimum of the objective function).
#' 
#' @details
#' EI selects the next design point as the one that maximizes the expected improvement:
#' \deqn{x_{next} = \arg\max_{x} \textrm{EI}(x) = \sigma(x)\left[u\Phi(u)+\phi(u)\right],}{x_next = argmax_x EI(x) = sigma(x)[u*Phi(u) + phi(u)],}
#' where \eqn{u = (y_{min}-\mu(x))/\sigma(x)}{u = (y_min - mu(x))/sigma(x)},
#' \eqn{\mu(x)}{mu(x)} and \eqn{\sigma(x)}{sigma(x)} are the kriging predictive mean and standard deviation,
#' and \eqn{\Phi}{Phi} and \eqn{\phi}{phi} are the standard normal CDF and PDF.
#' This is the most popular acquisition function for Bayesian optimization.
#'
#' The design space is specified via \code{design.parameters}, which must contain
#' \code{design.lower.bound} and \code{design.upper.bound} vectors.
#' During optimization, candidate points can be provided via \code{design.candidates}
#' in \code{\link{Get.Next.Design}}.
#' If the kriging object was fitted with \code{scale=TRUE}, design bounds,
#' candidates, new design points, predictions, returned designs, and returned data
#' are supplied and returned on the original input scale; scaling is handled
#' internally for the sequential criterion.
#'
#' @param design.parameters a list specifying the design space, containing:
#' \describe{
#'   \item{\code{design.lower.bound}}{a vector for the lower bound of the design space}
#'   \item{\code{design.upper.bound}}{a vector for the upper bound of the design space}
#' }
#' @param kriging (optional if kriging.parameters provided) a fitted kriging class object (e.g., from \link{Ordinary.Kriging} or \link{Fit.Kriging})
#' @param kriging.parameters (optional if kriging provided) a list of parameters to create kriging via \link{Fit.Kriging},
#'        must contain X, y and optionally other \link{Fit.Kriging} parameters
#' @param nlopt.parameters a list of parameters for NLopt optimization (used when \code{design.candidates} is not provided)
#' 
#' @return 
#' An EI Sequential Class Object.
#' 
#' @author 
#' Chaofan Huang, Difan Song, and V. Roshan Joseph
#' 
#' @references
#' Jones, D. R., Schonlau, M., & Welch, W. J. (1998). Efficient global optimization 
#' of expensive black-box functions. \emph{Journal of Global optimization}, 13(4), 455-492.
#' 
#' @seealso 
#' \link{PI.Sequential}, \link{UCB.Sequential}, \link{UDF.Sequential}, 
#' \link{Get.Next.Design}, \link{Add.Next.Design}, \link{Get.Best.Design}.
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
#' n <- 10
#' p <- 1
#' X <- matrix(runif(n), ncol=p)
#' y <- apply(X, 1, f)
#' newX <- matrix(seq(0, 1, length=1001), ncol=p)
#'
#' # create sequential
#' kriging <- Ordinary.Kriging(X, y, interpolation=TRUE, fit=TRUE,
#'                             kernel.parameters=list(type="Gaussian"))
#' design.parameters <- list(design.lower.bound=c(0), design.upper.bound=c(1))
#' sequential <- EI.Sequential(design.parameters, kriging)
#'
#' # sequential design loop with visualization
#' for (i in 1:5) {
#'   # visualize current kriging
#'   pred <- Predict.Sequential(sequential, newX)
#'   plot(newX, f(newX), "l", main=sprintf("EI Iteration %d", i))
#'   points(Get.Data(sequential)$X, Get.Data(sequential)$y, pch=16, col="blue")
#'   lines(newX, pred$mean, col="red", lty=2)
#'   lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#'   lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#'   # get next design and evaluate
#'   Xn <- Get.Next.Design(sequential)
#'   new.design <- list(X=Xn, y=f(Xn))
#'   points(new.design$X, new.design$y, pch=16, col="green")
#'   # add to kriging and update
#'   sequential <- Add.Next.Design(sequential, new.design)
#' }
#'
#' # final visualization
#' pred <- Predict.Sequential(sequential, newX)
#' plot(newX, f(newX), "l", main="EI Final")
#' points(Get.Data(sequential)$X, Get.Data(sequential)$y, pch=16, col="blue")
#' lines(newX, pred$mean, col="red", lty=2)
#' lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#' lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#' best.design <- Get.Best.Design(sequential)
#' points(best.design$design, best.design$objective, pch=16, col="green")
#' 
EI.Sequential <- function(design.parameters, kriging=NULL, kriging.parameters=NULL, nlopt.parameters=list()) {
  # check that at least one of kriging or kriging.parameters is provided
  if (is.null(kriging) & is.null(kriging.parameters))
    stop("Either kriging object or kriging.parameters must be provided.")
  
  # create kriging if kriging is not provided
  if (is.null(kriging)) {
    if (!is.list(kriging.parameters))
      stop("kriging.parameters must be a list.")
    if (!all(c("X", "y") %in% names(kriging.parameters)))
      stop("kriging.parameters must contain X and y.")
    kriging <- do.call(Fit.Kriging, kriging.parameters)
  }
  
  if (!grepl("^Rcpp.*Kriging$", class(kriging)))
    stop("An invalid kriging object is provided.")
  design.setup <- Get.Sequential.Design.Setup(design.parameters, kriging)
  sequential <- new(EISequential, kriging, design.setup$lower.bound, design.setup$upper.bound)
  sequential <- Set.Sequential.Scale.Info(sequential, design.setup$scale.info)
  sequential <- Set.Sequential.NLOPT.Parameters(sequential, nlopt.parameters)
  return (sequential)
}

#' @title 
#' PI Sequential Design
#' 
#' @description
#' This function creates a Probability of Improvement (PI) sequential design object 
#' for optimization (finding the minimum of the objective function).
#' 
#' @details
#' PI selects the next design point as the one that maximizes the probability of improvement:
#' \deqn{x_{next} = \arg\max_{x} \textrm{PI}(x) = \Phi(u),}{x_next = argmax_x PI(x) = Phi(u),}
#' where \eqn{u = (y_{min}-\mu(x))/\sigma(x)}{u = (y_min - mu(x))/sigma(x)},
#' \eqn{\mu(x)}{mu(x)} and \eqn{\sigma(x)}{sigma(x)} are the kriging predictive mean and standard deviation,
#' and \eqn{\Phi}{Phi} is the standard normal CDF.
#'
#' The design space is specified via \code{design.parameters}, which must contain
#' \code{design.lower.bound} and \code{design.upper.bound} vectors.
#' During optimization, candidate points can be provided via \code{design.candidates}
#' in \code{\link{Get.Next.Design}}.
#' If the kriging object was fitted with \code{scale=TRUE}, design bounds,
#' candidates, new design points, predictions, returned designs, and returned data
#' are supplied and returned on the original input scale; scaling is handled
#' internally for the sequential criterion.
#'
#' @param design.parameters a list specifying the design space, containing:
#' \describe{
#'   \item{\code{design.lower.bound}}{a vector for the lower bound of the design space}
#'   \item{\code{design.upper.bound}}{a vector for the upper bound of the design space}
#' }
#' @param kriging (optional if kriging.parameters provided) a fitted kriging class object (e.g., from \link{Ordinary.Kriging} or \link{Fit.Kriging})
#' @param kriging.parameters (optional if kriging provided) a list of parameters to create kriging via \link{Fit.Kriging},
#'        must contain X, y and optionally other \link{Fit.Kriging} parameters
#' @param nlopt.parameters a list of parameters for NLopt optimization (used when \code{design.candidates} is not provided)
#' 
#' @return 
#' A PI Sequential Class Object.
#' 
#' @author 
#' Chaofan Huang, Difan Song, and V. Roshan Joseph
#' 
#' @references
#' Kushner, H. J. (1964). A new method of locating the maximum point of an arbitrary multipeak 
#' curve in the presence of noise. \emph{Journal of Basic Engineering}, 86(1), 97-106.
#' 
#' @seealso 
#' \link{EI.Sequential}, \link{UCB.Sequential}, \link{UDF.Sequential}, 
#' \link{Get.Next.Design}, \link{Add.Next.Design}, \link{Get.Best.Design}.
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
#' n <- 10
#' p <- 1
#' X <- matrix(runif(n), ncol=p)
#' y <- apply(X, 1, f)
#' newX <- matrix(seq(0, 1, length=1001), ncol=p)
#'
#' # create sequential
#' kriging <- Ordinary.Kriging(X, y, interpolation=TRUE, fit=TRUE,
#'                             kernel.parameters=list(type="Gaussian"))
#' design.parameters <- list(design.lower.bound=c(0), design.upper.bound=c(1))
#' sequential <- PI.Sequential(design.parameters, kriging)
#'
#' # sequential design loop with visualization
#' for (i in 1:5) {
#'   # visualize current kriging
#'   pred <- Predict.Sequential(sequential, newX)
#'   plot(newX, f(newX), "l", main=sprintf("PI Iteration %d", i))
#'   points(Get.Data(sequential)$X, Get.Data(sequential)$y, pch=16, col="blue")
#'   lines(newX, pred$mean, col="red", lty=2)
#'   lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#'   lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#'   # get next design and evaluate
#'   Xn <- Get.Next.Design(sequential)
#'   new.design <- list(X=Xn, y=f(Xn))
#'   points(new.design$X, new.design$y, pch=16, col="green")
#'   # add to kriging and update
#'   sequential <- Add.Next.Design(sequential, new.design)
#' }
#'
#' # final visualization and best design
#' pred <- Predict.Sequential(sequential, newX)
#' plot(newX, f(newX), "l", main="PI Final")
#' points(Get.Data(sequential)$X, Get.Data(sequential)$y, pch=16, col="blue")
#' lines(newX, pred$mean, col="red", lty=2)
#' lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#' lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#' best.design <- Get.Best.Design(sequential)
#' points(best.design$design, best.design$objective, pch=16, col="green")
#' 
PI.Sequential <- function(design.parameters, kriging=NULL, kriging.parameters=NULL, nlopt.parameters=list()) {
  # check that at least one of kriging or kriging.parameters is provided
  if (is.null(kriging) & is.null(kriging.parameters))
    stop("Either kriging object or kriging.parameters must be provided.")
  
  # create kriging if kriging is not provided
  if (is.null(kriging)) {
    if (!is.list(kriging.parameters))
      stop("kriging.parameters must be a list.")
    if (!all(c("X", "y") %in% names(kriging.parameters)))
      stop("kriging.parameters must contain X and y.")
    kriging <- do.call(Fit.Kriging, kriging.parameters)
  }
  
  if (!grepl("^Rcpp.*Kriging$", class(kriging)))
    stop("An invalid kriging object is provided.")
  design.setup <- Get.Sequential.Design.Setup(design.parameters, kriging)
  sequential <- new(PISequential, kriging, design.setup$lower.bound, design.setup$upper.bound)
  sequential <- Set.Sequential.Scale.Info(sequential, design.setup$scale.info)
  sequential <- Set.Sequential.NLOPT.Parameters(sequential, nlopt.parameters)
  return (sequential)
}

#' @title 
#' UCB Sequential Design
#' 
#' @description
#' This function creates an Upper Confidence Bound (UCB) sequential design object 
#' for optimization (finding the minimum of the objective function).
#' 
#' @details
#' UCB selects the next design point as the one that minimizes the lower confidence bound:
#' \deqn{x_{next} = \arg\min_{x} \mu(x) - \lambda\sigma(x),}{x_next = argmin_x mu(x) - lambda*sigma(x),}
#' where \eqn{\mu(x)}{mu(x)} and \eqn{\sigma(x)}{sigma(x)} are the kriging predictive mean and standard deviation,
#' and \eqn{\lambda}{lambda} controls the exploration-exploitation tradeoff.
#' Larger \eqn{\lambda}{lambda} encourages more exploration.
#'
#' The design space is specified via \code{design.parameters}, which must contain
#' \code{design.lower.bound} and \code{design.upper.bound} vectors.
#' During optimization, candidate points can be provided via \code{design.candidates}
#' in \code{\link{Get.Next.Design}}.
#' If the kriging object was fitted with \code{scale=TRUE}, design bounds,
#' candidates, new design points, predictions, returned designs, and returned data
#' are supplied and returned on the original input scale; scaling is handled
#' internally for the sequential criterion.
#'
#' @param design.parameters a list specifying the design space, containing:
#' \describe{
#'   \item{\code{design.lower.bound}}{a vector for the lower bound of the design space}
#'   \item{\code{design.upper.bound}}{a vector for the upper bound of the design space}
#' }
#' @param kriging (optional if kriging.parameters provided) a fitted kriging class object (e.g., from \link{Ordinary.Kriging} or \link{Fit.Kriging})
#' @param kriging.parameters (optional if kriging provided) a list of parameters to create kriging via \link{Fit.Kriging},
#'        must contain X, y and optionally other \link{Fit.Kriging} parameters
#' @param lambda a positive scalar controlling the exploration-exploitation tradeoff (default: 2)
#' @param nlopt.parameters a list of parameters for NLopt optimization (used when \code{design.candidates} is not provided)
#' 
#' @return 
#' A UCB Sequential Class Object.
#' 
#' @author 
#' Chaofan Huang, Difan Song, and V. Roshan Joseph
#' 
#' @references
#' Srinivas, N., Krause, A., Kakade, S. M., & Seeger, M. (2009). Gaussian process optimization 
#' in the bandit setting: No regret and experimental design. \emph{arXiv preprint arXiv:0912.3995}.
#' 
#' @seealso 
#' \link{EI.Sequential}, \link{PI.Sequential}, \link{UDF.Sequential}, 
#' \link{Get.Next.Design}, \link{Add.Next.Design}, \link{Get.Best.Design}.
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
#' n <- 10
#' p <- 1
#' X <- matrix(runif(n), ncol=p)
#' y <- apply(X, 1, f)
#' newX <- matrix(seq(0, 1, length=1001), ncol=p)
#'
#' # create sequential
#' kriging <- Ordinary.Kriging(X, y, interpolation=TRUE, fit=TRUE,
#'                             kernel.parameters=list(type="Gaussian"))
#' design.parameters <- list(design.lower.bound=c(0), design.upper.bound=c(1))
#' sequential <- UCB.Sequential(design.parameters, kriging, lambda=2)
#'
#' # sequential design loop with visualization
#' for (i in 1:5) {
#'   # visualize current kriging
#'   pred <- Predict.Sequential(sequential, newX)
#'   plot(newX, f(newX), "l", main=sprintf("UCB Iteration %d", i))
#'   points(Get.Data(sequential)$X, Get.Data(sequential)$y, pch=16, col="blue")
#'   lines(newX, pred$mean, col="red", lty=2)
#'   lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#'   lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#'   # get next design and evaluate
#'   Xn <- Get.Next.Design(sequential)
#'   new.design <- list(X=Xn, y=f(Xn))
#'   points(new.design$X, new.design$y, pch=16, col="green")
#'   # add to kriging and update
#'   sequential <- Add.Next.Design(sequential, new.design)
#' }
#'
#' # final visualization and best design
#' pred <- Predict.Sequential(sequential, newX)
#' plot(newX, f(newX), "l", main="UCB Final")
#' points(Get.Data(sequential)$X, Get.Data(sequential)$y, pch=16, col="blue")
#' lines(newX, pred$mean, col="red", lty=2)
#' lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#' lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#' best.design <- Get.Best.Design(sequential)
#' points(best.design$design, best.design$objective, pch=16, col="green")
#' 
UCB.Sequential <- function(design.parameters, kriging=NULL, kriging.parameters=NULL, lambda=2, nlopt.parameters=list()) {
  # check that at least one of kriging or kriging.parameters is provided
  if (is.null(kriging) & is.null(kriging.parameters))
    stop("Either kriging object or kriging.parameters must be provided.")
  
  # create kriging if kriging is not provided
  if (is.null(kriging)) {
    if (!is.list(kriging.parameters))
      stop("kriging.parameters must be a list.")
    if (!all(c("X", "y") %in% names(kriging.parameters)))
      stop("kriging.parameters must contain X and y.")
    kriging <- do.call(Fit.Kriging, kriging.parameters)
  }
  
  if (!grepl("^Rcpp.*Kriging$", class(kriging)))
    stop("An invalid kriging object is provided.")
  design.setup <- Get.Sequential.Design.Setup(design.parameters, kriging)
  if (!(is.numeric(lambda) & lambda > 0))
    stop("lambda must be a positive number.")
  sequential <- new(UCBSequential, kriging, design.setup$lower.bound, design.setup$upper.bound, lambda)
  sequential <- Set.Sequential.Scale.Info(sequential, design.setup$scale.info)
  sequential <- Set.Sequential.NLOPT.Parameters(sequential, nlopt.parameters)
  return (sequential)
}

#' @title 
#' UDF Sequential Design
#' 
#' @description
#' This function creates a User Defined Function (UDF) sequential design object 
#' with a custom acquisition function. 
#' 
#' @details
#' The UDF sequential design allows users to specify a custom criterion function
#' for selecting the next design point. The criterion function receives the
#' predictive mean, predictive standard deviation, the current training data,
#' and any additional user-defined parameters.
#' If the kriging object was fitted with \code{scale=TRUE}, the criterion
#' function receives the current training \code{X} on the original input scale.
#'
#' The design space is specified via \code{design.parameters}, which must contain
#' \code{design.lower.bound} and \code{design.upper.bound} vectors.
#' During optimization, candidate points can be provided via \code{design.candidates}
#' in \code{\link{Get.Next.Design}}.
#' If the kriging object was fitted with \code{scale=TRUE}, design bounds,
#' candidates, new design points, predictions, returned designs, and returned data
#' are supplied and returned on the original input scale; scaling is handled
#' internally for the sequential criterion.
#'
#' @param design.parameters a list specifying the design space, containing:
#' \describe{
#'   \item{\code{design.lower.bound}}{a vector for the lower bound of the design space}
#'   \item{\code{design.upper.bound}}{a vector for the upper bound of the design space}
#' }
#' @param kriging (optional if kriging.parameters provided) a fitted kriging class object (e.g., from \link{Ordinary.Kriging} or \link{Fit.Kriging})
#' @param kriging.parameters (optional if kriging provided) a list of parameters to create kriging via \link{Fit.Kriging},
#'        must contain X, y and optionally other \link{Fit.Kriging} parameters
#' @param criterion.function an R function that takes (mean, sd, X, y, criterion.parameters) and returns a scalar criterion value to be maximized
#' @param criterion.parameters a list of additional parameters to pass to the criterion function
#' @param nlopt.parameters a list of parameters for NLopt optimization (used when \code{design.candidates} is not provided)
#' 
#' @return 
#' A UDF Sequential Class Object.
#' 
#' @author 
#' Chaofan Huang, Difan Song, and V. Roshan Joseph
#' 
#' @seealso 
#' \link{ALM.Sequential}, \link{EI.Sequential}, \link{PI.Sequential}, \link{UCB.Sequential}, 
#' \link{Predict.Sequential}, \link{Get.Next.Design}, \link{Add.Next.Design}.
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
#' n <- 10
#' p <- 1
#' X <- matrix(runif(n), ncol=p)
#' y <- apply(X, 1, f)
#' newX <- matrix(seq(0, 1, length=1001), ncol=p)
#'
#' # create sequential
#' kriging <- Ordinary.Kriging(X, y, interpolation=TRUE, fit=TRUE,
#'                             kernel.parameters=list(type="Gaussian"))
#' design.parameters <- list(design.lower.bound=c(0), design.upper.bound=c(1))
#'
#' # define custom criterion function (e.g., ALM-like using predictive SD)
#' my.criterion <- function(mean, sd, X, y, clist) { return(sd) }
#' sequential <- UDF.Sequential(design.parameters, kriging,
#'                              criterion.function=my.criterion,
#'                              criterion.parameters=list())
#'
#' # sequential design loop with visualization
#' for (i in 1:5) {
#'   # visualize current kriging
#'   pred <- Predict.Sequential(sequential, newX)
#'   plot(newX, f(newX), "l", main=sprintf("UDF Iteration %d", i))
#'   points(Get.Data(sequential)$X, Get.Data(sequential)$y, pch=16, col="blue")
#'   lines(newX, pred$mean, col="red", lty=2)
#'   lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#'   lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#'   # get next design and evaluate
#'   Xn <- Get.Next.Design(sequential)
#'   new.design <- list(X=Xn, y=f(Xn))
#'   points(new.design$X, new.design$y, pch=16, col="green")
#'   # add to kriging and update
#'   sequential <- Add.Next.Design(sequential, new.design)
#' }
#'
#' # final visualization
#' pred <- Predict.Sequential(sequential, newX)
#' plot(newX, f(newX), "l", main="UDF Final")
#' points(Get.Data(sequential)$X, Get.Data(sequential)$y, pch=16, col="blue")
#' lines(newX, pred$mean, col="red", lty=2)
#' lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#' lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#' 
UDF.Sequential <- function(design.parameters, kriging=NULL, kriging.parameters=NULL, criterion.function,
                            criterion.parameters=list(), nlopt.parameters=list()) {
  # check that at least one of kriging or kriging.parameters is provided
  if (is.null(kriging) & is.null(kriging.parameters))
    stop("Either kriging object or kriging.parameters must be provided.")
  
  # create kriging if kriging is not provided
  if (is.null(kriging)) {
    if (!is.list(kriging.parameters))
      stop("kriging.parameters must be a list.")
    if (!all(c("X", "y") %in% names(kriging.parameters)))
      stop("kriging.parameters must contain X and y.")
    kriging <- do.call(Fit.Kriging, kriging.parameters)
  }
  
  if (!grepl("^Rcpp.*Kriging$", class(kriging)))
    stop("An invalid kriging object is provided.")
  design.setup <- Get.Sequential.Design.Setup(design.parameters, kriging)
  if (!is.function(criterion.function))
    stop("criterion.function must be an R function.")
  if (!is.list(criterion.parameters))
    stop("criterion.parameters must be a list.")
  criterion.function <- Wrap.Scaled.Sequential.Criterion(criterion.function, design.setup$scale.info)
  sequential <- new(UDFSequential, kriging, design.setup$lower.bound, design.setup$upper.bound, criterion.function, criterion.parameters)
  sequential <- Set.Sequential.Scale.Info(sequential, design.setup$scale.info)
  sequential <- Set.Sequential.NLOPT.Parameters(sequential, nlopt.parameters)
  return (sequential)
}

#' @title
#' Add Next Design to Sequential Object
#'
#' @description
#' This function adds a new design point and its corresponding output to the sequential 
#' design object and updates kriging (optionally refit the kriging hyperparameters).
#'
#' @details
#' The new design is provided as a list containing:
#' \describe{
#'   \item{\code{X}}{a matrix or vector of new design point(s); if a vector, it will be converted to a 1×p matrix}
#'   \item{\code{y}}{a vector of corresponding output value(s)}
#' }
#' 
#' This function is typically used in a sequential design loop together with 
#' \code{\link{Get.Next.Design}}, allowing fine-grained control over each iteration.
#' If the sequential object was created from a scaled kriging object, \code{new.design$X}
#' should be provided on the original input scale and is scaled internally before
#' updating the kriging model.
#'
#' @param sequential a sequential design object (e.g., from \link{ALM.Sequential}, \link{EI.Sequential}, etc.)
#' @param new.design a list containing:
#'   \describe{
#'     \item{\code{X}}{matrix or vector of new design point(s) to add}
#'     \item{\code{y}}{vector of corresponding output value(s)}
#'   }
#' @param refit.kriging a logical indicating whether to refit the kriging hyperparameters
#'
#' @return
#' The sequential object (modified in-place).
#'
#' @author
#' Chaofan Huang, Difan Song, and V. Roshan Joseph
#'
#' @seealso
#' \link{Get.Next.Design}, \link{Get.Sequential}.
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
#' n <- 10
#' p <- 1
#' X <- matrix(runif(n), ncol=p)
#' y <- apply(X, 1, f)
#' newX <- matrix(seq(0, 1, length=1001), ncol=p)
#'
#' # create sequential
#' kriging <- Ordinary.Kriging(X, y, interpolation=TRUE, fit=TRUE,
#'                             kernel.parameters=list(type="Gaussian"))
#' design.parameters <- list(design.lower.bound=c(0), design.upper.bound=c(1))
#' sequential <- EI.Sequential(design.parameters, kriging)
#'
#' # sequential design loop with visualization
#' for (i in 1:5) {
#'   # visualize current kriging
#'   pred <- Predict.Sequential(sequential, newX)
#'   plot(newX, f(newX), "l", main=sprintf("EI Iteration %d", i))
#'   points(Get.Data(sequential)$X, Get.Data(sequential)$y, pch=16, col="blue")
#'   lines(newX, pred$mean, col="red", lty=2)
#'   lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#'   lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#'   # get next design and evaluate
#'   Xn <- Get.Next.Design(sequential)
#'   new.design <- list(X=Xn, y=f(Xn))
#'   points(new.design$X, new.design$y, pch=16, col="green")
#'   # add to kriging and update
#'   sequential <- Add.Next.Design(sequential, new.design)
#' }
#'
#' # final visualization
#' pred <- Predict.Sequential(sequential, newX)
#' plot(newX, f(newX), "l", main="EI Final")
#' points(Get.Data(sequential)$X, Get.Data(sequential)$y, pch=16, col="blue")
#' lines(newX, pred$mean, col="red", lty=2)
#' lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#' lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#' best.design <- Get.Best.Design(sequential)
#' points(best.design$design, best.design$objective, pch=16, col="green")
#'
Add.Next.Design <- function(sequential, new.design, refit.kriging=TRUE) {
  if (!grepl("^Rcpp.*Sequential$", class(sequential)))
    stop("An invalid sequential object is provided.")
  if (!is.list(new.design))
    stop("new.design must be a list containing X and y.")
  if (!all(c("X", "y") %in% names(new.design)))
    stop("new.design must contain both X and y components.")
  
  # Extract X and y from the list
  Xn <- new.design$X
  yn <- new.design$y
  
  # Ensure Xn is a matrix
  if (!is.matrix(Xn)) Xn <- matrix(Xn, nrow=1)
  # Ensure yn is a vector
  if (!is.vector(yn)) yn <- c(yn)
  # Ensure size and dimension are correct  
  p <- sequential$get_design_dimension()
  if (ncol(Xn) != p)
    stop(sprintf("Design dimension (%d) does not match sequential dimension (%d).", 
                 ncol(Xn), p))
  if (nrow(Xn) != length(yn))
    stop("Number of design points and outputs must match.")
  
  # Add data and refit kriging
  Xn <- Scale.Sequential.Input(sequential, Xn)
  sequential$add_data(Xn, yn)
  sequential$fit_kriging(refit.kriging)
  
  return(sequential)
}

#' @title
#' Get Next Design from Sequential
#'
#' @description
#' This function selects the next design point using the sequential criterion.
#'
#' @details
#' If \code{design.candidates} is provided, the next design is selected from the candidate set
#' by evaluating the criterion at each candidate point.
#' Otherwise, NLopt optimization is used within the design space bounds
#' specified in \code{design.parameters}.
#' If the sequential object was created from a scaled kriging object,
#' \code{design.candidates} should be provided on the original input scale and
#' the returned next design is also on the original input scale.
#'
#' @param sequential a sequential design object (e.g., from \link{ALM.Sequential}, \link{EI.Sequential}, etc.)
#' @param design.candidates (optional) a matrix of candidate design points; if provided,
#'        the next design is selected from this set instead of NLopt optimization
#'
#' @return
#' A vector representing the next design point.
#'
#' @author
#' Chaofan Huang, Difan Song, and V. Roshan Joseph
#'
#' @seealso
#' \link{Get.Sequential}, \link{Add.Next.Design}.
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
#' n <- 10
#' p <- 1
#' X <- matrix(runif(n), ncol=p)
#' y <- apply(X, 1, f)
#' newX <- matrix(seq(0, 1, length=1001), ncol=p)
#'
#' # fit kriging
#' kriging <- Ordinary.Kriging(X, y, interpolation=TRUE, fit=TRUE,
#'                             kernel.parameters=list(type="Gaussian"))
#' design.parameters <- list(design.lower.bound=c(0), design.upper.bound=c(1))
#' sequential <- EI.Sequential(design.parameters, kriging)
#'
#' # sequential design loop with visualization
#' for (i in 1:5) {
#'   # visualize current kriging
#'   pred <- Predict.Sequential(sequential, newX)
#'   plot(newX, f(newX), "l", main=sprintf("EI Iteration %d", i))
#'   points(Get.Data(sequential)$X, Get.Data(sequential)$y, pch=16, col="blue")
#'   lines(newX, pred$mean, col="red", lty=2)
#'   lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#'   lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#'   # get next design and evaluate
#'   Xn <- Get.Next.Design(sequential)
#'   new.design <- list(X=Xn, y=f(Xn))
#'   points(new.design$X, new.design$y, pch=16, col="green")
#'   # add to kriging and update
#'   sequential <- Add.Next.Design(sequential, new.design)
#' }
#'
#' # final visualization
#' pred <- Predict.Sequential(sequential, newX)
#' plot(newX, f(newX), "l", main="EI Final")
#' points(Get.Data(sequential)$X, Get.Data(sequential)$y, pch=16, col="blue")
#' lines(newX, pred$mean, col="red", lty=2)
#' lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#' lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#' best.design <- Get.Best.Design(sequential)
#' points(best.design$design, best.design$objective, pch=16, col="green")
#'
#' # alternatively can pass a set of candidate points for finding the next design
#' design.candidates <- matrix(seq(0,1,length=1001), ncol=p)
#' Get.Next.Design(sequential, design.candidates=design.candidates)
#'
Get.Next.Design <- function(sequential, design.candidates=NULL) {
  if (!grepl("^Rcpp.*Sequential$", class(sequential)))
    stop("An invalid sequential object is provided.")
  p <- sequential$get_design_dimension()
  design.lower.bound <- sequential$get_design_lower_bound()
  design.upper.bound <- sequential$get_design_upper_bound()
  original.design.lower.bound <- c(Unscale.Sequential.Input(sequential, matrix(design.lower.bound, nrow=1)))
  original.design.upper.bound <- c(Unscale.Sequential.Input(sequential, matrix(design.upper.bound, nrow=1)))
  if (!is.null(design.candidates)) {
    if (!is.matrix(design.candidates)) design.candidates <- as.matrix(design.candidates)
    if (ncol(design.candidates) != p)
      stop(sprintf("Dimension of design.candidates (%d) does not match with the sequential dimension (%d).",
                   ncol(design.candidates), p))
    if (any(sweep(design.candidates, 2L, original.design.lower.bound, "<")) ||
        any(sweep(design.candidates, 2L, original.design.upper.bound, ">")))
      stop("Design candidates must be within the design bounds.")
    design.candidates <- Scale.Sequential.Input(sequential, design.candidates)
    xn <- sequential$get_next_design_from_candidates(design.candidates)
  } else {
    x0 <- design.lower.bound + runif(p) * (design.upper.bound - design.lower.bound)
    xn <- sequential$get_next_design(x0)
  }
  xn <- matrix(xn, nrow=1, ncol=p)
  xn <- Unscale.Sequential.Input(sequential, xn)
  return (xn)
}

#' @title 
#' Predict from Sequential Design
#' 
#' @description
#' This function computes kriging predictions from the sequential design object at new input locations.
#'
#' @details
#' If the sequential object was created from a scaled kriging object, \code{newX}
#' should be provided on the original input scale and is scaled internally before
#' prediction.
#' 
#' @param sequential a sequential design object (e.g., from \link{ALM.Sequential}, \link{EI.Sequential}, etc.)
#' @param newX a matrix of new input locations for prediction
#' 
#' @return 
#' A list with components:
#' \describe{
#'   \item{\code{mean}}{the predictive mean}
#'   \item{\code{sd}}{the predictive standard deviation}
#' }
#' 
#' @author 
#' Chaofan Huang, Difan Song, and V. Roshan Joseph
#' 
#' @seealso
#' \link{Get.Sequential}.
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
#' n <- 10
#' p <- 1
#' X <- matrix(runif(n), ncol=p)
#' y <- apply(X, 1, f)
#' newX <- matrix(seq(0, 1, length=1001), ncol=p)
#'
#' # create sequential
#' kriging <- Ordinary.Kriging(X, y, interpolation=TRUE, fit=TRUE,
#'                             kernel.parameters=list(type="Gaussian"))
#' design.parameters <- list(design.lower.bound=c(0), design.upper.bound=c(1))
#' sequential <- EI.Sequential(design.parameters, kriging)
#' # get prediction and visualization
#' pred <- Predict.Sequential(sequential, newX)
#' plot(newX, f(newX), "l", main="EI Final")
#' points(Get.Data(sequential)$X, Get.Data(sequential)$y, pch=16, col="blue")
#' lines(newX, pred$mean, col="red", lty=2)
#' lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#' lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#' 
Predict.Sequential <- function(sequential, newX) {
  if (!grepl("^Rcpp.*Sequential$", class(sequential)))
    stop("An invalid sequential object is provided.")
  if (!is.matrix(newX)) newX <- as.matrix(newX)
  if (ncol(newX) != sequential$get_design_dimension())
    stop(sprintf("Dimension of newX (%d) does not match with the sequential dimension (%d).", 
                 ncol(newX), sequential$get_design_dimension()))
  newX <- Scale.Sequential.Input(sequential, newX)
  return (sequential$predict(newX))
}

#' @title 
#' Get Best Design from Optimization Sequential
#' 
#' @description
#' This function returns the best (minimum) objective value and its corresponding design 
#' from an optimization sequential design object (e.g., \link{EI.Sequential}, \link{PI.Sequential}, \link{UCB.Sequential}).
#'
#' @details
#' If the sequential object was created from a scaled kriging object, the returned
#' best design is converted back to the original input scale.
#' 
#' @param sequential an optimization sequential design object
#' 
#' @return 
#' A list with components:
#' \describe{
#'   \item{\code{design}}{the design (input) corresponding to the best objective}
#'   \item{\code{objective}}{the best (minimum) objective value}
#' }
#' 
#' @author 
#' Chaofan Huang, Difan Song, and V. Roshan Joseph
#' 
#' @seealso 
#' \link{EI.Sequential}, \link{PI.Sequential}, \link{UCB.Sequential}.
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
#' n <- 10
#' p <- 1
#' X <- matrix(runif(n), ncol=p)
#' y <- apply(X, 1, f)
#' newX <- matrix(seq(0, 1, length=1001), ncol=p)
#'
#' # fit kriging
#' kriging <- Ordinary.Kriging(X, y, interpolation=TRUE, fit=TRUE,
#'                             kernel.parameters=list(type="Gaussian"))
#' design.parameters <- list(design.lower.bound=c(0), design.upper.bound=c(1))
#' sequential <- EI.Sequential(design.parameters, kriging)
#'
#' # sequential design loop with visualization
#' for (i in 1:5) {
#'   # visualize current kriging
#'   pred <- Predict.Sequential(sequential, newX)
#'   plot(newX, f(newX), "l", main=sprintf("EI Iteration %d", i))
#'   points(Get.Data(sequential)$X, Get.Data(sequential)$y, pch=16, col="blue")
#'   lines(newX, pred$mean, col="red", lty=2)
#'   lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#'   lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#'   # get next design and evaluate
#'   Xn <- Get.Next.Design(sequential)
#'   new.design <- list(X=Xn, y=f(Xn))
#'   points(new.design$X, new.design$y, pch=16, col="green")
#'   # add to kriging and update
#'   sequential <- Add.Next.Design(sequential, new.design)
#' }
#'
#' # final visualization
#' pred <- Predict.Sequential(sequential, newX)
#' plot(newX, f(newX), "l", main="EI Final")
#' points(Get.Data(sequential)$X, Get.Data(sequential)$y, pch=16, col="blue")
#' lines(newX, pred$mean, col="red", lty=2)
#' lines(newX, pred$mean-2*pred$sd, col="red", lty=3)
#' lines(newX, pred$mean+2*pred$sd, col="red", lty=3)
#' best.design <- Get.Best.Design(sequential)
#' points(best.design$design, best.design$objective, pch=16, col="green")
#' 
Get.Best.Design <- function(sequential) {
  if (!grepl("^Rcpp.*Sequential$", class(sequential)))
    stop("An invalid sequential object is provided.")
  if (!inherits(sequential, "Rcpp_OptimizationSequential") & 
      !inherits(sequential, "Rcpp_EISequential") & 
      !inherits(sequential, "Rcpp_PISequential") & 
      !inherits(sequential, "Rcpp_UCBSequential"))
    stop("Get.Best.Design is only available for optimization sequential objects (EI, PI, UCB).")
  best.design <- matrix(sequential$get_best_design(), nrow=1)
  best.design <- c(Unscale.Sequential.Input(sequential, best.design))
  return (list(
    design = best.design,
    objective = sequential$get_best_objective()
  ))
}

#' @title
#' Get Data from Sequential
#'
#' @description
#' This function retrieves the current training data (X and y) from a sequential design object.
#'
#' @details
#' The sequential design object maintains the current training data internally.
#' This function provides access to both the design matrix X and the corresponding output vector y.
#' If the sequential object was created from a scaled kriging object, \code{X} is
#' converted back to the original input scale.
#'
#' @param sequential a sequential design object (e.g., from \link{ALM.Sequential}, \link{EI.Sequential}, etc.)
#'
#' @return
#' A list with components:
#' \describe{
#'   \item{\code{X}}{the design matrix (training inputs)}
#'   \item{\code{y}}{the output vector (training responses)}
#' }
#'
#' @author
#' Chaofan Huang, Difan Song, and V. Roshan Joseph
#'
#' @seealso
#' \link{Get.Sequential}.
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
#' n <- 10
#' p <- 1
#' X <- matrix(runif(n), ncol=p)
#' y <- apply(X, 1, f)
#' newX <- matrix(seq(0, 1, length=1001), ncol=p)
#'
#' # create sequential
#' kriging <- Ordinary.Kriging(X, y, interpolation=TRUE, fit=TRUE,
#'                             kernel.parameters=list(type="Gaussian"))
#' design.parameters <- list(design.lower.bound=c(0), design.upper.bound=c(1))
#' sequential <- EI.Sequential(design.parameters, kriging)
#' Get.Data(sequential)
#'
Get.Data <- function(sequential) {
  if (!grepl("^Rcpp.*Sequential$", class(sequential)))
    stop("An invalid sequential object is provided.")
  data <- sequential$get_data()
  data$X <- Unscale.Sequential.Input(sequential, data$X)
  return (data)
}
