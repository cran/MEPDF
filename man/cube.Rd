\name{cube}
\alias{cube}
\title{cube}
\usage{
cube(data,mx,mn,grid.sizes)
}
\description{
Single grid size empirical density function. To be used to construct the epdf function.
}
\arguments{
\item{data}{N-dimensional data set.}
\item{mx}{Upper cropping point of the data.}
\item{mn}{Lower cropping point of the data.}
\item{grid.sizes}{Vector of grid sizes for the grid.}
}
\examples{
library("pracma")
library("plyr")

data<-cbind(rnorm(1000),rnorm(1000))

pdf<-cube(data,mx=c(1,1),mn=c(-1,-1),grid.sizes = c(0.01,0.01))
}
