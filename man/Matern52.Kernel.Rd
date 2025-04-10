% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernel.R
\name{Matern52.Kernel}
\alias{Matern52.Kernel}
\title{Matern(5/2) Kernel}
\usage{
Matern52.Kernel(lengthscale)
}
\arguments{
\item{lengthscale}{a vector for the positive length scale parameters}
}
\value{
A Matern(5/2) Kernel Class Object.
}
\description{
This function specifies the Matern kernel with smoothness parameter \eqn{\nu}{nu}=5/2.
}
\details{
The Matern(5/2) kernel is given by 
\deqn{k(r)=(1+\sqrt{5}r+5r^2/3)\exp(-\sqrt{5}r),}{k(r) = (1+sqrt(3)r+5r^2/3)exp(-sqrt(5)r),}
where 
\deqn{r(x,x^{\prime})=\sqrt{\sum_{i=1}^{p}\left(\frac{x_{i}-x_{i}^{\prime}}{l_{i}}\right)^2}}{r(x,x') = sqrt(sum_{i=1}^{p} [(x_i - x'_i) / l_i]^2)} 
is the euclidean distance between \eqn{x}{x} and \eqn{x^{\prime}}{x'} weighted by
the length scale parameters \eqn{l_{i}}{l_i}'s.
}
\examples{
n <- 5
p <- 3
X <- matrix(rnorm(n*p), ncol=p)
lengthscale <- c(1:p)

# approach 1
kernel <- Matern52.Kernel(lengthscale)
Evaluate.Kernel(kernel, X)

# approach 2
kernel <- Get.Kernel(lengthscale, type="Matern52")
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
