\name{normkernel}
\alias{normkernel}
\title{normkernel}
\usage{
normkernel(x,H)
}
\description{
Kernel function based on the Epanechnikov kernel.
}
\arguments{
\item{x}{Evaluation point.}
\item{H}{Positive-definite, symmetric matrix as bandwidth.}
}
\examples{
normkernel(c(1,1),H = diag(2))
}
