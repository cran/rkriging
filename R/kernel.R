#' Kernel Class
#'
#' @useDynLib Kernel, .registration = TRUE
#' @import Rcpp, RcppEigen, BH, nloptr
#' 
#' @noRd
Rcpp::loadModule("Kernel", TRUE)

#' @title 
#' Kernel
#' 
#' @description
#' This function provides a common interface to specify various kernels.
#' See arguments section for the available kernels in this package. 
#' 
#' @param lengthscale a vector for the positive length scale parameters
#' @param type kernel type: Gaussian, RQ, Matern12, Matern32, Matern52, Matern, UDF, MultiplicativeRQ, MultiplicativeMatern, MultiplicativeUDF 
#' @param parameters a list of parameters required for the specific kernel
#' 
#' @return 
#' A Kernel Class Object.
#' 
#' @author 
#' Chaofan Huang and V. Roshan Joseph
#' 
#' @references
#' Duvenaud, D. (2014). \emph{The kernel cookbook: Advice on covariance functions}.
#' 
#' Rasmussen, C. E. & Williams, C. K. (2006). \emph{Gaussian Processes for Machine Learning}. The MIT Press.
#' 
#' @seealso 
#' \link{Evaluate.Kernel}, \link{Gaussian.Kernel}, \link{RQ.Kernel}, 
#' \link{Matern12.Kernel}, \link{Matern32.Kernel}, \link{Matern52.Kernel}
#' \link{Matern.Kernel}, \link{UDF.Kernel},
#' \link{MultiplicativeRQ.Kernel}, \link{MultiplicativeMatern.Kernel}, \link{MultiplicativeUDF.Kernel}.
#' 
#' @export
#' 
#' @examples 
#' n <- 5
#' p <- 3
#' X <- matrix(rnorm(n*p), ncol=p)
#' lengthscale <- c(1:p)
#' 
#' # Gaussian 
#' kernel <- Get.Kernel(lengthscale, type="Gaussian")
#' Evaluate.Kernel(kernel, X)
#' 
#' # Rational Quadratic (RQ)
#' kernel <- Get.Kernel(lengthscale, type="RQ", parameters=list(alpha=1))
#' Evaluate.Kernel(kernel, X) 
#' 
#' # Matern(1/2)
#' kernel <- Get.Kernel(lengthscale, type="Matern12")
#' Evaluate.Kernel(kernel, X) 
#' 
#' # Matern(3/2)
#' kernel <- Get.Kernel(lengthscale, type="Matern32")
#' Evaluate.Kernel(kernel, X) 
#' 
#' # Matern(5/2)
#' kernel <- Get.Kernel(lengthscale, type="Matern52")
#' Evaluate.Kernel(kernel, X) 
#' 
#' # Generalized Matern
#' kernel <- Get.Kernel(lengthscale, type="Matern", parameters=list(nu=2.01))
#' Evaluate.Kernel(kernel, X) 
#' 
#' # User Defined Function (UDF) Kernel
#' kernel.function <- function(sqdist) {return (exp(-sqrt(sqdist)))} 
#' kernel <- Get.Kernel(lengthscale, type="UDF", 
#'                      parameters=list(kernel.function=kernel.function))
#' Evaluate.Kernel(kernel, X) 
#' 
#' # Multiplicative Rational Quadratic (RQ)
#' kernel <- Get.Kernel(lengthscale, type="MultiplicativeRQ", parameters=list(alpha=1))
#' Evaluate.Kernel(kernel, X) 
#' 
#' # Multiplicative Generalized Matern
#' kernel <- Get.Kernel(lengthscale, type="MultiplicativeMatern", parameters=list(nu=2.01))
#' Evaluate.Kernel(kernel, X)
#' 
#' # Multiplicative User Defined Function (UDF)
#' kernel.function <- function(sqdist) {return (exp(-sqrt(sqdist)))} 
#' kernel <- Get.Kernel(lengthscale, type="MultiplicativeUDF", 
#'                      parameters=list(kernel.function=kernel.function))
#' Evaluate.Kernel(kernel, X) 
#' 
Get.Kernel <- function(lengthscale, type, parameters=list()) {
  if (!is.character(type))
      stop("kernel type must be a string.")
  if (!(is.vector(lengthscale) & is.numeric(lengthscale) & all(lengthscale>0)))
    stop("kernel lengthscale must be a vector of positive number.")
  if (type == "Gaussian") {
    kernel <- new(GaussianKernel, lengthscale)
  } else if (type == "RQ") {
    if (!("alpha" %in% names(parameters)))
      stop("alpha must be provided for RQ Kernel.")
    alpha <- parameters[["alpha"]]
    if (!(is.numeric(alpha) & alpha>0))
      stop("alpha must be positive for RQ Kernel.")
    kernel <- new(RQKernel, lengthscale, alpha)
  } else if (type == "Matern12") {
    kernel <- new(Matern12Kernel, lengthscale)
  } else if (type == "Matern32") {
    kernel <- new(Matern32Kernel, lengthscale)
  } else if (type == "Matern52") {
    kernel <- new(Matern52Kernel, lengthscale)
  } else if (type == "Matern") {
    if (!("nu" %in% names(parameters)))
      stop("nu must be provided for Matern Kernel.")
    nu <- parameters[["nu"]]
    if (!(is.numeric(nu) & nu>0))
      stop("nu must be positive for Matern Kernel.")
    kernel <- new(MaternKernel, lengthscale, nu)
  } else if (type == "UDF") {
    if (!("kernel.function" %in% names(parameters)))
      stop("kernel.function must be provided for UDF Kernel.")
    kernel.function <- parameters[["kernel.function"]]
    if (!is.function(kernel.function))
      stop("kernel.function must be a function for UDF Kernel.")
    kernel <- new(UDFKernel, lengthscale, kernel.function)
  } else if (type == "MultiplicativeRQ") {
    if (!("alpha" %in% names(parameters)))
      stop("alpha must be provided for MultiplicativeRQ Kernel.")
    alpha <- parameters[["alpha"]]
    if (!(is.numeric(alpha) & alpha>0))
      stop("alpha must be positive for MultiplicativeRQ Kernel.")
    kernel <- new(MultiplicativeRQKernel, lengthscale, alpha)
  } else if (type == "MultiplicativeMatern") {
    if (!("nu" %in% names(parameters)))
      stop("nu must be provided for MultiplicativeMatern Kernel.")
    nu <- parameters[["nu"]]
    if (!(is.numeric(nu) & nu>0))
      stop("nu must be positive for MultiplicativeMatern Kernel.")
    kernel <- new(MultiplicativeMaternKernel, lengthscale, nu)
  } else if (type == "MultiplicativeUDF") {
    if (!("kernel.function" %in% names(parameters)))
      stop("kernel.function must be provided for MultiplicativeUDF Kernel.")
    kernel.function <- parameters[["kernel.function"]]
    if (!is.function(kernel.function))
      stop("kernel.function must be a function for MultiplicativeUDF Kernel.")
    kernel <- new(MultiplicativeUDFKernel, lengthscale, kernel.function)
  } else {
    stop(sprintf("No Kernel type %s.", type))
  }
  return (kernel)
}

#' @title 
#' Gaussian Kernel
#' 
#' @description
#' This function specifies the Gaussian / Squared Exponential (SE) / Radial Basis Function (RBF) kernel.
#' 
#' @details
#' The Gaussian kernel is given by 
#' \deqn{k(r)=\exp(-r^2/2),}{k(r) = exp(-r^2/2),}
#' where 
#' \deqn{r(x,x^{\prime})=\sqrt{\sum_{i=1}^{p}\left(\frac{x_{i}-x_{i}^{\prime}}{l_{i}}\right)^2}}{r(x,x') = sqrt(sum_{i=1}^{p} [(x_i - x'_i) / l_i]^2)} 
#' is the euclidean distance between \eqn{x}{x} and \eqn{x^{\prime}}{x'} weighted by
#' the length scale parameters \eqn{l_{i}}{l_i}'s.
#' 
#' @param lengthscale a vector for the positive length scale parameters
#' 
#' @return 
#' A Gaussian Kernel Class Object.
#' 
#' @author 
#' Chaofan Huang and V. Roshan Joseph
#' 
#' @references
#' Duvenaud, D. (2014). \emph{The kernel cookbook: Advice on covariance functions}.
#' 
#' Rasmussen, C. E. & Williams, C. K. (2006). \emph{Gaussian Processes for Machine Learning}. The MIT Press.
#' 
#' @seealso 
#' \link{Get.Kernel}, \link{Evaluate.Kernel}.
#' 
#' @export
#' 
#' @examples 
#' n <- 5
#' p <- 3
#' X <- matrix(rnorm(n*p), ncol=p)
#' lengthscale <- c(1:p)
#' 
#' # approach 1
#' kernel <- Gaussian.Kernel(lengthscale)
#' Evaluate.Kernel(kernel, X)
#' 
#' # approach 2
#' kernel <- Get.Kernel(lengthscale, type="Gaussian")
#' Evaluate.Kernel(kernel, X) 
#' 
Gaussian.Kernel <- function(lengthscale) {
  if (!(is.vector(lengthscale) & is.numeric(lengthscale) & all(lengthscale>0)))
    stop("lengthscale must be a vector of positive number.")
  kernel <- new(GaussianKernel, lengthscale)
  return (kernel)
}

#' @title 
#' Rational Quadratic (RQ) Kernel
#' 
#' @description
#' This function specifies the Rational Quadratic (RQ) kernel.
#' 
#' @details
#' The Rational Quadratic (RQ) kernel is given by 
#' \deqn{k(r;\alpha)=\left(1+\frac{r^2}{2\alpha}\right)^{-\alpha},}{k(r;alpha) = (1+[r^2/(2*alpha)])^{-alpha},}
#' where \eqn{\alpha}{alpha} is the scale mixture parameter and 
#' \deqn{r(x,x^{\prime})=\sqrt{\sum_{i=1}^{p}\left(\frac{x_{i}-x_{i}^{\prime}}{l_{i}}\right)^2}}{r(x,x') = sqrt(sum_{i=1}^{p} [(x_i - x'_i) / l_i]^2)} 
#' is the euclidean distance between \eqn{x}{x} and \eqn{x^{\prime}}{x'} weighted by
#' the length scale parameters \eqn{l_{i}}{l_i}'s.
#' As \eqn{\alpha\to\infty}{alpha goes to infinity}, it converges to the \link{Gaussian.Kernel}.
#' 
#' @param lengthscale a vector for the positive length scale parameters
#' @param alpha a positive scalar for the scale mixture parameter that controls the relative weighting of large-scale and small-scale variations
#' 
#' @return 
#' A Rational Quadratic (RQ) Kernel Class Object.
#' 
#' @author 
#' Chaofan Huang and V. Roshan Joseph
#' 
#' @references
#' Duvenaud, D. (2014). \emph{The kernel cookbook: Advice on covariance functions}.
#' 
#' Rasmussen, C. E. & Williams, C. K. (2006). \emph{Gaussian Processes for Machine Learning}. The MIT Press.
#' 
#' @seealso 
#' \link{MultiplicativeRQ.Kernel}, \link{Get.Kernel}, \link{Evaluate.Kernel}.
#' 
#' @export
#' 
#' @examples 
#' n <- 5
#' p <- 3
#' X <- matrix(rnorm(n*p), ncol=p)
#' lengthscale <- c(1:p)
#' 
#' # approach 1
#' kernel <- RQ.Kernel(lengthscale, alpha=1)
#' Evaluate.Kernel(kernel, X)
#' 
#' # approach 2
#' kernel <- Get.Kernel(lengthscale, type="RQ", parameters=list(alpha=1))
#' Evaluate.Kernel(kernel, X) 
#' 
RQ.Kernel <- function(lengthscale, alpha=1) {
  if (!(is.vector(lengthscale) & is.numeric(lengthscale) & all(lengthscale>0)))
    stop("lengthscale must be a vector of positive number.")
  if (!(is.numeric(alpha) & alpha>0))
    stop("alpha must be positive for RQ Kernel.")
  kernel <- new(RQKernel, lengthscale, alpha)
  return (kernel)
}

#' @title 
#' Matern(1/2) Kernel
#' 
#' @description
#' This function specifies the Matern kernel with smoothness parameter \eqn{\nu}{nu}=1/2. 
#' It is also known as the Exponential kernel. 
#' 
#' @details
#' The Matern(1/2) kernel is given by 
#' \deqn{k(r)=\exp(-r),}{k(r) = exp(-r),}
#' where 
#' \deqn{r(x,x^{\prime})=\sqrt{\sum_{i=1}^{p}\left(\frac{x_{i}-x_{i}^{\prime}}{l_{i}}\right)^2}}{r(x,x') = sqrt(sum_{i=1}^{p} [(x_i - x'_i) / l_i]^2)} 
#' is the euclidean distance between \eqn{x}{x} and \eqn{x^{\prime}}{x'} weighted by
#' the length scale parameters \eqn{l_{i}}{l_i}'s.
#' 
#' @param lengthscale a vector for the positive length scale parameters
#' 
#' @return 
#' A Matern(1/2) Kernel Class Object.
#' 
#' @author 
#' Chaofan Huang and V. Roshan Joseph
#' 
#' @references
#' Duvenaud, D. (2014). \emph{The kernel cookbook: Advice on covariance functions}.
#' 
#' Rasmussen, C. E. & Williams, C. K. (2006). \emph{Gaussian Processes for Machine Learning}. The MIT Press.
#' 
#' @seealso 
#' \link{Matern.Kernel}, \link{Get.Kernel}, \link{Evaluate.Kernel}.
#' 
#' @export
#' 
#' @examples 
#' n <- 5
#' p <- 3
#' X <- matrix(rnorm(n*p), ncol=p)
#' lengthscale <- c(1:p)
#' 
#' # approach 1
#' kernel <- Matern12.Kernel(lengthscale)
#' Evaluate.Kernel(kernel, X)
#' 
#' # approach 2
#' kernel <- Get.Kernel(lengthscale, type="Matern12")
#' Evaluate.Kernel(kernel, X) 
#' 
Matern12.Kernel <- function(lengthscale) {
  if (!(is.vector(lengthscale) & is.numeric(lengthscale) & all(lengthscale>0)))
    stop("lengthscale must be a vector of positive number.")
  kernel <- new(Matern12Kernel, lengthscale)
  return (kernel)
}

#' @title 
#' Matern(3/2) Kernel
#' 
#' @description
#' This function specifies the Matern kernel with smoothness parameter \eqn{\nu}{nu}=3/2. 
#' 
#' @details
#' The Matern(3/2) kernel is given by 
#' \deqn{k(r)=(1+\sqrt{3}r)\exp(-\sqrt{3}r),}{k(r) = (1+sqrt(3)r)exp(-sqrt(3)r),}
#' where 
#' \deqn{r(x,x^{\prime})=\sqrt{\sum_{i=1}^{p}\left(\frac{x_{i}-x_{i}^{\prime}}{l_{i}}\right)^2}}{r(x,x') = sqrt(sum_{i=1}^{p} [(x_i - x'_i) / l_i]^2)} 
#' is the euclidean distance between \eqn{x}{x} and \eqn{x^{\prime}}{x'} weighted by
#' the length scale parameters \eqn{l_{i}}{l_i}'s.
#' 
#' @param lengthscale a vector for the positive length scale parameters
#' 
#' @return 
#' A Matern(3/2) Kernel Class Object.
#' 
#' @author 
#' Chaofan Huang and V. Roshan Joseph
#' 
#' @references
#' Duvenaud, D. (2014). \emph{The kernel cookbook: Advice on covariance functions}.
#' 
#' Rasmussen, C. E. & Williams, C. K. (2006). \emph{Gaussian Processes for Machine Learning}. The MIT Press.
#' 
#' @seealso 
#' \link{Matern.Kernel}, \link{Get.Kernel}, \link{Evaluate.Kernel}.
#' 
#' @export
#' 
#' @examples 
#' n <- 5
#' p <- 3
#' X <- matrix(rnorm(n*p), ncol=p)
#' lengthscale <- c(1:p)
#' 
#' # approach 1
#' kernel <- Matern32.Kernel(lengthscale)
#' Evaluate.Kernel(kernel, X)
#' 
#' # approach 2
#' kernel <- Get.Kernel(lengthscale, type="Matern32")
#' Evaluate.Kernel(kernel, X) 
#' 
Matern32.Kernel <- function(lengthscale) {
  if (!(is.vector(lengthscale) & is.numeric(lengthscale) & all(lengthscale>0)))
    stop("lengthscale must be a vector of positive number.")
  kernel <- new(Matern32Kernel, lengthscale)
  return (kernel)
}

#' @title 
#' Matern(5/2) Kernel
#' 
#' @description
#' This function specifies the Matern kernel with smoothness parameter \eqn{\nu}{nu}=5/2. 
#' 
#' @details
#' The Matern(5/2) kernel is given by 
#' \deqn{k(r)=(1+\sqrt{5}r+5r^2/3)\exp(-\sqrt{5}r),}{k(r) = (1+sqrt(3)r+5r^2/3)exp(-sqrt(5)r),}
#' where 
#' \deqn{r(x,x^{\prime})=\sqrt{\sum_{i=1}^{p}\left(\frac{x_{i}-x_{i}^{\prime}}{l_{i}}\right)^2}}{r(x,x') = sqrt(sum_{i=1}^{p} [(x_i - x'_i) / l_i]^2)} 
#' is the euclidean distance between \eqn{x}{x} and \eqn{x^{\prime}}{x'} weighted by
#' the length scale parameters \eqn{l_{i}}{l_i}'s.
#' 
#' @param lengthscale a vector for the positive length scale parameters
#' 
#' @return 
#' A Matern(5/2) Kernel Class Object.
#' 
#' @author 
#' Chaofan Huang and V. Roshan Joseph
#' 
#' @references
#' Duvenaud, D. (2014). \emph{The kernel cookbook: Advice on covariance functions}.
#' 
#' Rasmussen, C. E. & Williams, C. K. (2006). \emph{Gaussian Processes for Machine Learning}. The MIT Press.
#' 
#' @seealso 
#' \link{Matern.Kernel}, \link{Get.Kernel}, \link{Evaluate.Kernel}.
#' 
#' @export
#' 
#' @examples 
#' n <- 5
#' p <- 3
#' X <- matrix(rnorm(n*p), ncol=p)
#' lengthscale <- c(1:p)
#' 
#' # approach 1
#' kernel <- Matern52.Kernel(lengthscale)
#' Evaluate.Kernel(kernel, X)
#' 
#' # approach 2
#' kernel <- Get.Kernel(lengthscale, type="Matern52")
#' Evaluate.Kernel(kernel, X) 
#' 
Matern52.Kernel <- function(lengthscale) {
  if (!(is.vector(lengthscale) & is.numeric(lengthscale) & all(lengthscale>0)))
    stop("lengthscale must be a vector of positive number.")
  kernel <- new(Matern52Kernel, lengthscale)
  return (kernel)
}

#' @title 
#' Generalized Matern Kernel
#' 
#' @description
#' This function specifies the (Generalized) Matern kernel with any smoothness parameter \eqn{\nu}{nu}.
#' 
#' @details
#' The Generalized Matern kernel is given by 
#' \deqn{k(r;\nu)=\frac{2^{1-\nu}}{\Gamma(\nu)}(\sqrt{2\nu}r)^{\nu}K_{\nu}(\sqrt{2\nu}r),}{k(r;nu) = 2^(1-nu) / Gamma(nu) * (sqrt(2*nu)*r)^(nu) * K_nu(sqrt(2*nu)*r),}
#' where \eqn{\nu}{nu} is the smoothness parameter, 
#' \eqn{K_{\nu}}{K_nu} is the modified Bessel function, 
#' \eqn{\Gamma}{Gamma} is the gamma function, 
#' and \deqn{r(x,x^{\prime})=\sqrt{\sum_{i=1}^{p}\left(\frac{x_{i}-x_{i}^{\prime}}{l_{i}}\right)^2}}{r(x,x') = sqrt(sum_{i=1}^{p} [(x_i - x'_i) / l_i]^2)} 
#' is the euclidean distance between \eqn{x}{x} and \eqn{x^{\prime}}{x'} weighted by
#' the length scale parameters \eqn{l_{i}}{l_i}'s. 
#' As \eqn{\nu\to\infty}{nu goes to infinity}, it converges to the \link{Gaussian.Kernel}. 
#' 
#' @param lengthscale a vector for the positive length scale parameters
#' @param nu a positive scalar parameter that controls the smoothness
#' 
#' @return 
#' A Generalized Matern Kernel Class Object.
#' 
#' @author 
#' Chaofan Huang and V. Roshan Joseph
#' 
#' @references
#' Duvenaud, D. (2014). \emph{The kernel cookbook: Advice on covariance functions}.
#' 
#' Rasmussen, C. E. & Williams, C. K. (2006). \emph{Gaussian Processes for Machine Learning}. The MIT Press.
#' 
#' @seealso 
#' \link{Matern12.Kernel}, \link{Matern32.Kernel}, \link{Matern52.Kernel}, 
#' \link{MultiplicativeMatern.Kernel}, \link{Get.Kernel}, \link{Evaluate.Kernel}.
#' 
#' @export
#' 
#' @examples 
#' n <- 5
#' p <- 3
#' X <- matrix(rnorm(n*p), ncol=p)
#' lengthscale <- c(1:p)
#' 
#' # approach 1
#' kernel <- Matern.Kernel(lengthscale, nu=2.01)
#' Evaluate.Kernel(kernel, X)
#' 
#' # approach 2
#' kernel <- Get.Kernel(lengthscale, type="Matern", parameters=list(nu=2.01))
#' Evaluate.Kernel(kernel, X) 
#' 
Matern.Kernel <- function(lengthscale, nu=2.01) {
  if (!(is.vector(lengthscale) & is.numeric(lengthscale) & all(lengthscale>0)))
    stop("lengthscale must be a vector of positive number.")
  if (!(is.numeric(nu) & nu>0))
    stop("nu must be positive for Matern Kernel.")
  kernel <- new(MaternKernel, lengthscale, nu)
  return (kernel)
}

#' @title 
#' User Defined Function (UDF) Kernel
#' 
#' @description
#' This function specifies a kernel with the user defined R function. 
#' 
#' @details
#' The User Defined Function (UDF) kernel is given by
#' \deqn{k(r) = f(r)}{k(r) = f(r)}
#' where \eqn{f}{f} is the user defined kernel function that takes \eqn{r^2}{r^2} as input, 
#' where \deqn{r(x,x^{\prime})=\sqrt{\sum_{i=1}^{p}\left(\frac{x_{i}-x_{i}^{\prime}}{l_{i}}\right)^2},}{r(x,x') = sqrt(sum_{i=1}^{p} [(x_i - x'_i) / l_i]^2),} 
#' is the euclidean distance between \eqn{x}{x} and \eqn{x^{\prime}}{x'} weighted by
#' the length scale parameters \eqn{l_{i}}{l_i}'s.
#' 
#' @param lengthscale a vector for the positive length scale parameters
#' @param kernel.function user defined kernel function
#' 
#' @return 
#' A User Defined Function (UDF) Kernel Class Object.
#' 
#' @references
#' Duvenaud, D. (2014). \emph{The kernel cookbook: Advice on covariance functions}.
#' 
#' Rasmussen, C. E. & Williams, C. K. (2006). \emph{Gaussian Processes for Machine Learning}. The MIT Press.
#' 
#' @author 
#' Chaofan Huang and V. Roshan Joseph
#' 
#' @seealso 
#' \link{MultiplicativeUDF.Kernel}, \link{Get.Kernel}, \link{Evaluate.Kernel}.
#' 
#' @export
#' 
#' @examples 
#' n <- 5
#' p <- 3
#' X <- matrix(rnorm(n*p), ncol=p)
#' lengthscale <- c(1:p)
#' 
#' kernel.function <- function(sqdist) {return (exp(-sqrt(sqdist)))} 
#' 
#' # approach 1
#' kernel <- UDF.Kernel(lengthscale, kernel.function=kernel.function)
#' Evaluate.Kernel(kernel, X)
#' 
#' # approach 2
#' kernel <- Get.Kernel(lengthscale, type="UDF", 
#'                      parameters=list(kernel.function=kernel.function))
#' Evaluate.Kernel(kernel, X) 
#' 
UDF.Kernel <- function(lengthscale, kernel.function) {
  if (!is.function(kernel.function))
    stop("kernel.function must be a function for UDF Kernel.")
  kernel <- new(UDFKernel, lengthscale, kernel.function)
  return (kernel)
}

#' @title 
#' Multiplicative Rational Quadratic (RQ) Kernel
#' 
#' @description
#' This function specifies the Multiplicative Rational Quadratic (RQ) kernel.
#' 
#' @details
#' The Multiplicative Rational Quadratic (RQ) kernel is given by 
#' \deqn{k(r;\alpha)=\prod_{i=1}^{p}\left(1+\frac{r_{i}^2}{2\alpha}\right)^{-\alpha},}{k(r; alpha) = prod_{i=1}^{p}(1+[r_i^2/(2*alpha)])^{-alpha},}
#' where \eqn{\alpha}{alpha} is the scale mixture parameter and 
#' \deqn{r_{i}(x,x^{\prime})=\sqrt{\left(\frac{x_{i}-x_{i}^{\prime}}{l_{i}}\right)^2}}{r_i(x,x') = sqrt([(x_i - x'_i) / l_i]^2)} 
#' is the dimension-wise euclidean distances between \eqn{x}{x} and \eqn{x^{\prime}}{x'} weighted by
#' the length scale parameters \eqn{l_{i}}{l_i}'s.
#' 
#' @param lengthscale a vector for the positive length scale parameters
#' @param alpha a positive scalar for the scale mixture parameter that controls the relative weighting of large-scale and small-scale variations
#' 
#' @return 
#' A Multiplicative Rational Quadratic (RQ) Kernel Class Object.
#' 
#' @author 
#' Chaofan Huang and V. Roshan Joseph
#' 
#' @references
#' Duvenaud, D. (2014). \emph{The kernel cookbook: Advice on covariance functions}.
#' 
#' Rasmussen, C. E. & Williams, C. K. (2006). \emph{Gaussian Processes for Machine Learning}. The MIT Press.
#' 
#' @seealso 
#' \link{RQ.Kernel}, \link{Get.Kernel}, \link{Evaluate.Kernel}.
#' 
#' @export
#' 
#' @examples 
#' n <- 5
#' p <- 3
#' X <- matrix(rnorm(n*p), ncol=p)
#' lengthscale <- c(1:p)
#' 
#' # approach 1
#' kernel <- MultiplicativeRQ.Kernel(lengthscale, alpha=1)
#' Evaluate.Kernel(kernel, X)
#' 
#' # approach 2
#' kernel <- Get.Kernel(lengthscale, type="MultiplicativeRQ", parameters=list(alpha=1))
#' Evaluate.Kernel(kernel, X) 
#' 
MultiplicativeRQ.Kernel <- function(lengthscale, alpha=1) {
  if (!(is.vector(lengthscale) & is.numeric(lengthscale) & all(lengthscale>0)))
    stop("lengthscale must be a vector of positive number.")
  if (!(is.numeric(alpha) & alpha>0))
    stop("alpha must be positive for MultiplicativeRQ Kernel.")
  kernel <- new(MultiplicativeRQKernel, lengthscale, alpha)
  return (kernel)
}

#' @title 
#' Multiplicative Generalized Matern Kernel
#' 
#' @description
#' This function specifies the Multiplicative Generalized Matern kernel.
#' 
#' @details
#' The Multiplicative Generalized Matern kernel is given by 
#' \deqn{k(r;\nu)=\prod_{i=1}^{p}\frac{2^{1-\nu}}{\Gamma(\nu)}(\sqrt{2\nu}r_{i})^{\nu}K_{\nu}(\sqrt{2\nu}r_{i}),}{k(r; nu) = prod_{i=1}^{p} 2^(1-nu) / Gamma(nu) * (sqrt(2*nu)*r_i)^(nu) * K_nu(sqrt(2*nu)*r_i),}
#' where \eqn{\nu}{nu} is the smoothness parameter, 
#' \eqn{K_{\nu}}{K_nu} is the modified Bessel function, 
#' \eqn{\Gamma}{Gamma} is the gamma function, 
#' and \deqn{r_{i}(x,x^{\prime})=\sqrt{\left(\frac{x_{i}-x_{i}^{\prime}}{l_{i}}\right)^2}}{r_i(x,x') = sqrt([(x_i - x'_i) / l_i]^2)} 
#' is the dimension-wise euclidean distances between \eqn{x}{x} and \eqn{x^{\prime}}{x'} weighted by
#' the length scale parameters \eqn{l_{i}}{l_i}'s. 
#' 
#' @param lengthscale a vector for the positive length scale parameters
#' @param nu a positive scalar parameter that controls the smoothness
#' 
#' @return 
#' A Multiplicative Generalized Matern Kernel Class Object.
#' 
#' @author 
#' Chaofan Huang and V. Roshan Joseph
#' 
#' @references
#' Duvenaud, D. (2014). \emph{The kernel cookbook: Advice on covariance functions}.
#' 
#' Rasmussen, C. E. & Williams, C. K. (2006). \emph{Gaussian Processes for Machine Learning}. The MIT Press.
#' 
#' @seealso 
#' \link{Matern.Kernel}, \link{Get.Kernel}, \link{Evaluate.Kernel}.
#' 
#' @export
#' 
#' @examples 
#' n <- 5
#' p <- 3
#' X <- matrix(rnorm(n*p), ncol=p)
#' lengthscale <- c(1:p)
#' 
#' # approach 1
#' kernel <- MultiplicativeMatern.Kernel(lengthscale, nu=2.01)
#' Evaluate.Kernel(kernel, X)
#' 
#' # approach 2
#' kernel <- Get.Kernel(lengthscale, type="MultiplicativeMatern", parameters=list(nu=2.01))
#' Evaluate.Kernel(kernel, X) 
#' 
MultiplicativeMatern.Kernel <- function(lengthscale, nu=2.01) {
  if (!(is.vector(lengthscale) & is.numeric(lengthscale) & all(lengthscale>0)))
    stop("lengthscale must be a vector of positive number.")
  if (!(is.numeric(nu) & nu>0))
    stop("nu must be positive for Matern Kernel.")
  kernel <- new(MultiplicativeMaternKernel, lengthscale, nu)
  return (kernel)
}

#' @title 
#' Multiplicative User Defined Function (UDF) Kernel
#' 
#' @description
#' This function specifies the Multiplicative kernel with the user defined R function. 
#' 
#' @details
#' The Multiplicative User Defined Function (UDF) kernel is given by 
#' \deqn{k(r)=\prod_{i=1}^{p}f(r_{i}),}{k(r) = prod_{i=1}^{p} f(r_i),}
#' where \eqn{f}{f} is the user defined kernel function that takes \eqn{r_{i}^2}{r_i^2} as input, 
#' where \deqn{r_{i}(x,x^{\prime})=\sqrt{\left(\frac{x_{i}-x_{i}^{\prime}}{l_{i}}\right)^2}}{r_i(x,x') = sqrt([(x_i - x'_i) / l_i]^2)} 
#' is the dimension-wise euclidean distances between \eqn{x}{x} and \eqn{x^{\prime}}{x'} weighted by
#' the length scale parameters \eqn{l_{i}}{l_i}'s. 
#' 
#' @param lengthscale a vector for the positive length scale parameters
#' @param kernel.function user defined kernel function
#' 
#' @return 
#' A Multiplicative User Defined Function (UDF) Kernel Class Object.
#' 
#' @references
#' Duvenaud, D. (2014). \emph{The kernel cookbook: Advice on covariance functions}.
#' 
#' Rasmussen, C. E. & Williams, C. K. (2006). \emph{Gaussian Processes for Machine Learning}. The MIT Press.
#' 
#' @author 
#' Chaofan Huang and V. Roshan Joseph
#' 
#' @seealso 
#' \link{UDF.Kernel}, \link{Get.Kernel}, \link{Evaluate.Kernel}.
#' 
#' @export
#' 
#' @examples 
#' n <- 5
#' p <- 3
#' X <- matrix(rnorm(n*p), ncol=p)
#' lengthscale <- c(1:p)
#' 
#' kernel.function <- function(sqdist) {return (exp(-sqrt(sqdist)))} 
#' 
#' # approach 1
#' kernel <- MultiplicativeUDF.Kernel(lengthscale, kernel.function=kernel.function)
#' Evaluate.Kernel(kernel, X)
#' 
#' # approach 2
#' kernel <- Get.Kernel(lengthscale, type="MultiplicativeUDF", 
#'                      parameters=list(kernel.function=kernel.function))
#' Evaluate.Kernel(kernel, X) 
#' 
MultiplicativeUDF.Kernel <- function(lengthscale, kernel.function) {
  if (!is.function(kernel.function))
    stop("kernel.function must be a function for the UDF Kernel.")
  kernel <- new(MultiplicativeUDFKernel, lengthscale, kernel.function)
  return (kernel)
}

#' @title 
#' Evaluate Kernel
#' 
#' @description
#' This function computes the kernel (correlation) matrix.
#' Given the kernel class object and the input data \eqn{X}{X} of size n, 
#' this function computes the corresponding \eqn{n\times n}{n x n} kernel (correlation) matrix. 
#' 
#' @param kernel a kernel class object
#' @param X input data
#' 
#' @return 
#' The kernel (correlation) matrix of X evaluated by the kernel.
#' 
#' @author 
#' Chaofan Huang and V. Roshan Joseph
#' 
#' @seealso 
#' \link{Get.Kernel}.
#' 
#' @export
#' 
#' @examples 
#' n <- 5
#' p <- 3
#' X <- matrix(rnorm(n*p), ncol=p)
#' lengthscale <- c(1:p)
#' 
#' kernel <- Gaussian.Kernel(lengthscale)
#' Evaluate.Kernel(kernel, X)
#' 
Evaluate.Kernel <- function(kernel, X) {
  if (!grepl("^Rcpp.*Kernel$", class(kernel)))
    stop("An invalid kernel object is provided.")
  if (!is.matrix(X)) X <- as.matrix(X)
  if (ncol(X) != kernel$get_dimension())
    stop(sprintf("Dimension of X (%d) does not match with the kernel dimension (%d).", 
                 ncol(X), kernel$get_dimension()))
  return (kernel$compute(X))
}
