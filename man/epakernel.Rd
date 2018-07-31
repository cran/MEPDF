\name{epakernel}
\alias{epakernel}
\title{epakernel}
\usage{
epakernel(x,H)
}
\description{
Kernel function based on the normal distribution.
}
\arguments{
\item{x}{Evaluation point.}
\item{H}{Positive-definite, symmetric matrix as bandwidth.}
}
\examples{
epakernel(c(1,1),H = diag(2))
}
