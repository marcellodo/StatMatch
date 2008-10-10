% --- Source file: create_fused.Rd ---
\name{create.fused}
\alias{create.fused}
\title{Creates a matched (synthetic) dataset}

\description{
Creates a \emph{synthetic} data frame after the statistical matching of two data sources at \emph{micro} level.
}

\usage{
create.fused(data.rec, data.don, mtc.ids, 
                z.vars, dup.x=FALSE, match.vars=NULL)  
}

\arguments{

\item{data.rec}{
A matrix or data frame that has been used as \emph{recipient} in the statistical matching application. 
}

\item{data.don}{
A matrix or data frame that has been used as the \emph{donor} in the statistical matching application.
}

\item{mtc.ids}{
A matrix with two columns. Each row must contain the name or the index of the recipient record (row) in \code{data.don} and the name or the index of the corresponding donor record (row) in \code{data.don}. Note that this type of matrix is returned by the functions \code{\link[StatMatch]{NND.hotdeck}} and \code{\link[StatMatch]{RANDwNND.hotdeck}}.
}

\item{z.vars}{
A character vector with the name of the variables available only in \code{data.don} that should be \dQuote{donated} to \code{data.rec}.
}

\item{dup.x}{
Logical. When \code{TRUE} the values of the matching variables in \code{data.don} are also \dQuote{donated} to \code{data.rec}. The name of the matching variables has to be specified with the argument \code{match.vars}. To avoid confusion, the matching variables added to \code{data.rec} are renamed by adding the suffix \dQuote{don}. By default \code{dup.x=FALSE}.
}

\item{match.vars}{
A character vector with the names of the matching variables. It has to be specified only when \code{dup.x=TRUE}.
}


}


  
\details{
This function allows to create the synthetic (or fused) data set after the application of a statistical matching in a  \emph{micro} framework. For details D'Orazio \emph{et al.} (2006).
} 


\value{

The \code{data.rec} data frame with the \code{z.vars} filled in and, when \code{dup.x=TRUE}, with the values of the matching variables for the donor records.

}
\references{
D'Orazio, M., Di Zio, M. and Scanu, M. (2006). \emph{Statistical Matching: Theory and Practice.} Wiley, Chichester.
}

\author{
 Marcello D'Orazio \email{madorazi@istat.it} 
}

\seealso{ 
\code{\link[StatMatch]{NNN.hotdeck}}
\code{\link[StatMatch]{RANDwNNN.hotdeck}}

}

\examples{

lab <- c(1:15, 51:65, 101:115)
iris.rec <- iris[lab, c(1:3,5)]  # recipient data.frame
iris.don <- iris[-lab, c(1:2,4:5)] # donor data.frame

# Now iris.rec and iris.don have the variables
# "Sepal.Length", "Sepal.Width"  and "Species"
# in common.
#  "Petal.Length" is available only in iris.rec
#  "Petal.Width"  is available only in iris.don

# find the closest donors using NND hot deck;
# distances are computed on "Sepal.Length" and "Sepal.Width"

out.NND <- NND.hotdeck(data.rec=iris.rec, data.don=iris.don,
            match.vars=c("Sepal.Length", "Sepal.Width"), don.class="Species")

# create synthetic data.set, without the duplication of the matching variables

fused.0 <- create.fused(data.rec=iris.rec, data.don=iris.don, 
            mtc.ids=out.NND$mtc.ids, z.vars="Petal.Width")

# create synthetic data.set, with the "duplication" of the matching variables

fused.1 <- create.fused(data.rec=iris.rec, data.don=iris.don,
            mtc.ids=out.NND$mtc.ids, z.vars="Petal.Width", dup.x=TRUE,
            match.vars=c("Sepal.Length", "Sepal.Width"))
}
\keyword{nonparametric}