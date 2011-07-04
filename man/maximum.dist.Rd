
\name{maximum.dist}
\alias{maximum.dist}
\title{Computes the Maximum  Distance}

\description{
  This function computes the Maximum distance (or \eqn{L^\infty}{L^Inf} norm) among units in a dataset or among observations in two distinct datasets.
}

\usage{
maximum.dist(data.x, data.y=data.x, rank=FALSE)
}

\arguments{

\item{data.x}{
  A matrix or a data frame containing variables that should be used in the computation of the distance.  Only continuous variables are allowed.  Missing values (\code{NA}) are not allowed.

When only \code{data.x} is supplied, the distances between rows of \code{data.x} is computed. 
}

\item{data.y}{
   A numeric matrix or data frame with the same variables, of the same type, as those in \code{data.x} (only continuous variables are allowed).  Dissimilarities between rows of \code{data.x} and rows of \code{data.y} will be computed.  If not provided, by default it is assumed \code{data.y=data.x} and only dissimilarities between rows of \code{data.x} will be computed. 
}

\item{rank}{
Logical, when \code{TRUE} the original values are substituted by their ranks divided by the number of values plus one (following suggestion in Kovar et al. 1988).  This rank transformation permits to remove the effect of different scales on the distance computation.  When computing ranks the tied observations assume the average of their position 
(\code{ties.method = "average"} in calling the \code{\link[base]{rank}} function). 

}
}

  
\details{
This function computes the \eqn{L^\infty}{L^Inf} distance also know as minimax distance.  In practice the distance among two records is the maximum of the absolute differences among the observed variables:

\deqn{d(i,j) = max \left( \left|x_{1i}-x_{1j} \right|, \left|x_{2i}-x_{2j} \right|,\ldots,\left|x_{Ki}-x_{Kj} \right| \right)}{d(i,j) =max(|x_1i-x_1j|,|x_2i-x_2j|,...,|x_Ki-x_Kj|) }

When \code{rank=TRUE} the original values are substituted by their ranks divided by the number of values plus one (following suggestion in Kovar et al. 1988). 

}

\value{
A \code{matrix} object with distances among rows of \code{data.x} and those of \code{data.y}.
}

\references{
Kovar, J.G., MacMillan, J. and Whitridge, P. (1988). \dQuote{Overview and strategy for the Generalized Edit and Imputation System}. Statistics Canada, Methodology Branch Working Paper No. BSMD 88-007 E/F.
}

\author{
 Marcello D'Orazio \email{madorazi@istat.it} 
}

\seealso{ 
\code{\link[base]{rank}}, 
}

\examples{

md1 <- maximum.dist(iris[1:10,1:4])
md2 <- maximum.dist(iris[1:10,1:4], rank=TRUE)

md3 <- maximum.dist(data.x=iris[1:50,1:4], data.y=iris[51:100,1:4])
md4 <- maximum.dist(data.x=iris[1:50,1:4], data.y=iris[51:100,1:4], rank=TRUE)

}

\keyword{multivariate}
