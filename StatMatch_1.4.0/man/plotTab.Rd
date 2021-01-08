% --- Source file: fact2dummy.Rd ---
\name{plotTab}
\alias{plotTab}
\title{graphical comparison of the estimated distributions for the same continuous variable.}

\description{
Compares graphically the estimated distributions for the same categorical variable(s) using data coming from two different data sources.
}

\usage{
plotTab(data.A, data.B, xlab.A, xlab.B=NULL, w.A=NULL, w.B=NULL)
}

\arguments{

\item{data.A}{
  	A dataframe or matrix containing the variable of interest \code{xlab.A} and eventual associated survey weights \code{w.A}. 
}

\item{data.B}{
  	A dataframe or matrix containing the variable of interest \code{xlab.B} and eventual associated survey weights \code{w.B}. 
}

\item{xlab.A}{
  	Character string providing the name(s) of one or more variables in \code{data.A} whose (joint) distribution should be represented graphically and compared with that estimated from \code{data.B}.
}

\item{xlab.B}{
  	Character string providing the name(s) of one or more variables in \code{data.A} whose (joint) distribution should be represented graphically and compared with that estimated from \code{data.A}. If \code{xlab.B=NULL} (default) then it assumed \code{xlab.B=xlab.A}.
}


\item{w.A}{
   Character string providing the name of the optional weighting variable in \code{data.A} that, in case, should be used to estimate the distribution of \code{xlab.A}
}

\item{w.B}{
   Character string providing the name of the optional weighting variable in \code{data.B} that, in case, should be used to estimate the distribution of \code{xlab.B}
}

}
  
\details{
This function compares graphically the (joint) distribution of the same variables but estimated from data coming from two different data sources. The graphical comparison is done using barcharts.
} 


\value{

The required graphical representation is drawn using the \pkg{ggplot2} facilities.

}

\author{
 Marcello D'Orazio \email{mdo.statmatch@gmail.com} 
}

\seealso{
\code{\link[StatMatch]{comp.prop}}
}

\examples{

# plotTab(data.A = samp.A, data.B = samp.B, xlab.A="edu7", w.A = "ww")
# plotTab(data.A = samp.A, data.B = samp.B, xlab.A=c("urb", "sex"), w.A = "ww", w.B="ww")

}

\keyword{multivariate}
\keyword{dplot}