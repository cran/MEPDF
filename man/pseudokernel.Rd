\name{pseudokernel}
\alias{pseudokernel}
\title{pseudokernel}
\usage{
pseudokernel(data,mn,mx,grid.sizes,rings)
}
\description{
Single grid size empirical density function. Offers a modified kernel density approach via the optional argument "rings"
}
\arguments{
\item{data}{N-dimensional data set.}
\item{mx}{Upper cropping point of the data.}
\item{mn}{Lower cropping point of the data.}
\item{grid.sizes}{Vector of grid sizes for the grid.}
\item{rings}{Number of additional cell grid smoothing layers.}
}
\examples{
data<-mvtnorm::rmvnorm(n = 100,mean = c(0,0),sigma = diag(2))

pdf<-pseudokernel(data = data,mn = c(-2,-2),mx = c(2,2),grid.sizes = c(0.1,0.1),rings = 1)
}
