\name{ekde}
\alias{ekde}
\title{ekde}
\usage{
ekde(x,data,H,rule,kernel)
}
\description{
Kernel function based on the normal distribution.
}
\arguments{
\item{x}{Evaluation point.}
\item{data}{Input data.}
\item{H}{Positive-definite, symmetric matrix as bandwidth.}
\item{rule}{In absence of a bandwidth matrix a rule-of-thumb can be chosen, either the "silverman" or "scott" rule.}
\item{kernel}{The kernel function of choice.}
}
\examples{
library("pracma")
library("plyr")

data<-cbind(rnorm(1000),rnorm(1000))
pdf<-ekde(x = 0,data = data, rule = "silverman",kernel = normkernel)
}
