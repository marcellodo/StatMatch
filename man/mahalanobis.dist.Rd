\name{mahalanobis.dist}
\alias{mahalanobis.dist}
\title{Computes the Mahalanobis Distance}

\description{
  This function computes the Mahalanobis distance among units in a dataset or among observations in two distinct datasets.
}

\usage{
mahalanobis.dist(data.x, data.y=NULL, vc=NULL)
}

\arguments{

\item{data.x}{
  A matrix or a data frame containing variables that should be used in the computation of the distance.  Only continuous variables are allowed.  Missing values (\code{NA}) are not allowed.

When only \code{data.x} is supplied, the distances between rows of \code{data.x} is computed. 
}

\item{data.y}{
   A numeric matrix or data frame with the same variables, of the same type, as those in \code{data.x} (only continuous variables are allowed).  Dissimilarities between rows of \code{data.x} and rows of \code{data.y} will be computed.  If not provided, by default it is assumed \code{data.y=data.x} and only dissimilarities between rows of \code{data.x} will be computed. 
}

\item{vc}{
Covariance matrix that should be used in distance computation.  If it is not supplied (\code{vc=NULL}) it is estimated from the input data.  In particular, when \code{vc=NULL} and only \code{data.x} is supplied then the covariance matrix is estimated from \code{data.x} (\code{var(data.x)}).  On the contrary when \code{vc=NULL} and both \code{data.x} and \code{data.y} are available then the covariance matrix is estimated on the joined data sets (i.e. \code{var(rbind(data.x, data.y))}).
}

}


\details{
This function computes the Mahalanobis distance:

\deqn{d(i,j)=\sqrt{(x_i - x_j)^T S^{-1} (x_i - x_j)}}{d(i,j) = ((x_i - x_j)^T * S^(-1) * (x_i - x_j) )^(1/2)}

When \code{vc=NULL} the covariance matrix S is estimated from the available data (see argument \code{vc} for details) otherwise the one supplied via the argument \code{vc} is used. 

}

\value{
A \code{matrix} object with distances among rows of \code{data.x} and those of \code{data.y}.
}

\references{
Mahalanobis, P C (1936) \dQuote{On the generalised distance in statistics}. Proceedings of the National Institute of Sciences of India 2, pp. 49-55.
}

\author{
 Marcello D'Orazio \email{madorazi@istat.it} 
}

\seealso{ 
\code{\link[stats]{mahalanobis}} 
}

\examples{

md1 <- mahalanobis.dist(iris[1:6,1:4])
md2 <- mahalanobis.dist(data.x=iris[1:6,1:4], data.y=iris[51:60, 1:4])

vv <- var(iris[,1:4])
md1a <- mahalanobis.dist(data.x=iris[1:6,1:4], vc=vv)
md2a <- mahalanobis.dist(data.x=iris[1:6,1:4], data.y=iris[51:60, 1:4], vc=vv)

}

\keyword{multivariate}
