% --- Source file: fact2dummy.Rd ---
\name{fact2dummy}
\alias{fact2dummy}
\title{Transforms a categorical variable in a set of dummy variables}

\description{
Transforms a factor or more factors contained in a data frame in a set of dummy variables.
}

\usage{
fact2dummy(data, all=TRUE, lab="x") 
}

\arguments{

\item{data}{
  	A factor or a data frame that contains one or more factors (columns whose class is \dQuote{factor} or \dQuote{ordered}) that have to be substituted by dummy variables. 
}

\item{all}{
   Logical. When \code{all=TRUE} (default) the output matrix will contain as many dummy variables as the number of the levels of the factor variable. On the contrary, when \code{all=FALSE}, the dummy variable related to the last level of the factor is dropped.
}

\item{lab}{
   A character string with the name of the variable to be pasted with its levels. This is used only when  \code{data} is a factor. By default it is set to \dQuote{x}.
}

}
  
\details{
This function substitutes categorical variables  in the input data frame (columns whose class is \dQuote{factor} or \dQuote{ordered}) with the corresponding dummy variables.    
} 


\value{

A matrix with the dummy variables instead of initial factor variables.

}

\author{
 Marcello D'Orazio \email{mdo.statmatch@gmail.com} 
}

\seealso{ 
\code{\link[StatMatch]{gower.dist}}
}

\examples{

x <- runif(5)
y <- factor(c(1,2,1,2,2))
z <- ordered(c(1,2,3,2,2))
xyz <- data.frame(x,y,z)
fact2dummy(xyz)

fact2dummy(xyz, all=FALSE)


#example with iris data frame
str(iris)
ir.mat <- fact2dummy(iris)
head(ir.mat)

}

\keyword{multivariate}
\keyword{cluster}