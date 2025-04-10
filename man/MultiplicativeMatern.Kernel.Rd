% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernel.R
\name{MultiplicativeMatern.Kernel}
\alias{MultiplicativeMatern.Kernel}
\title{Multiplicative Generalized Matern Kernel}
\usage{
MultiplicativeMatern.Kernel(lengthscale, nu = 2.01)
}
\arguments{
\item{lengthscale}{a vector for the positive length scale parameters}

\item{nu}{a positive scalar parameter that controls the smoothness}
}
\value{
A Multiplicative Generalized Matern Kernel Class Object.
}
\description{
This function specifies the Multiplicative Generalized Matern kernel.
}
\details{
The Multiplicative Generalized Matern kernel is given by 
\deqn{k(r;\nu)=\prod_{i=1}^{p}\frac{2^{1-\nu}}{\Gamma(\nu)}(\sqrt{2\nu}r_{i})^{\nu}K_{\nu}(\sqrt{2\nu}r_{i}),}{k(r; nu) = prod_{i=1}^{p} 2^(1-nu) / Gamma(nu) * (sqrt(2*nu)*r_i)^(nu) * K_nu(sqrt(2*nu)*r_i),}
where \eqn{\nu}{nu} is the smoothness parameter, 
\eqn{K_{\nu}}{K_nu} is the modified Bessel function, 
\eqn{\Gamma}{Gamma} is the gamma function, 
and \deqn{r_{i}(x,x^{\prime})=\sqrt{\left(\frac{x_{i}-x_{i}^{\prime}}{l_{i}}\right)^2}}{r_i(x,x') = sqrt([(x_i - x'_i) / l_i]^2)} 
is the dimension-wise euclidean distances between \eqn{x}{x} and \eqn{x^{\prime}}{x'} weighted by
the length scale parameters \eqn{l_{i}}{l_i}'s.
}
\examples{
n <- 5
p <- 3
X <- matrix(rnorm(n*p), ncol=p)
lengthscale <- c(1:p)

# approach 1
kernel <- MultiplicativeMatern.Kernel(lengthscale, nu=2.01)
Evaluate.Kernel(kernel, X)

# approach 2
kernel <- Get.Kernel(lengthscale, type="MultiplicativeMatern", parameters=list(nu=2.01))
Evaluate.Kernel(kernel, X) 

}
\references{
Duvenaud, D. (2014). \emph{The kernel cookbook: Advice on covariance functions}.

Rasmussen, C. E. & Williams, C. K. (2006). \emph{Gaussian Processes for Machine Learning}. The MIT Press.
}
\seealso{
\link{Matern.Kernel}, \link{Get.Kernel}, \link{Evaluate.Kernel}.
}
\author{
Chaofan Huang and V. Roshan Joseph
}
