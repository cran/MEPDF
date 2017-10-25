\name{epdf}
\alias{epdf}
\title{epdf}
\usage{
epdf(data,max.corner,min.corner,main.gridsize,rescubes)
}
\description{
Assembles main grid and superimposes grids of different resolution.
}
\arguments{
\item{data}{N-dimensional data set.}
\item{max.corner}{Upper cropping point of the data.}
\item{min.corner}{Lower cropping point of the data.}
\item{main.gridsize}{Vector of grid sizes for the main grid.}
\item{rescubes}{List of upper and lower cropping points, as well as grid sizes for cubes that are to be superimposed.}
}
\examples{
library("pracma")
library("plyr")

a<-list(c(-1,-1),c(1,1),c(0.01,0.01))
b<-list(c(-2,-2),c(2,2),c(0.02,0.02))
cubes<-list(a,b)

min.corner = c(-4,-4)
max.corner= c(4,4)
main.gridsize = c(0.05,0.05)

# Data & Density
data<-cbind(rnorm(1000),rnorm(1000))

pdf<-epdf(data,max.corner,min.corner,main.gridsize,rescubes = cubes)

}
