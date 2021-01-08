% --- Source file: fact2dummy.Rd ---
\name{plotBounds}
\alias{plotBounds}
\title{Graphical representation of the uncertainty bounds estimated through the \code{\link[StatMatch]{Frechet.bounds.cat}}  function}

\description{
The function uses the output of the function \code{\link[StatMatch]{Frechet.bounds.cat}} to produce a basic graphical representation of the uncertainty bounds related to the contingency table of Y vs. Z.
}

\usage{
plotBounds(outFB)
}

\arguments{

\item{outFB}{
  	the list provided in output from \code{\link[StatMatch]{Frechet.bounds.cat}}. 
}

}
  
\details{
This function represents graphically the uncertainty bounds estimated by the function \code{\link[StatMatch]{Frechet.bounds.cat}} for each relative frequency in the contingency table of Y vs. Z.
the dotted line indicates the width of the bounds estimated without conditioning on the Xs (the size is reported in parenthesis below the line). The full line indicates the width of the estimated bounds conditional on the Xs (expected conditional Frechet bounds for the relative frequencies in the contingency table of Y vs. Z (size reported below the line, not in the the parenthesis).
Not that when the X are not used it is drawn only the width of the unconditional bounds and the size is shown below the line.

The figure on the top od the line indicated the estimated relative frequency under the assumption of independence between Y and Z conditional on one or more X variables (Conditional Independence Assumption, CIA; for details see help pages of \code{\link[StatMatch]{Frechet.bounds.cat}}), otherwise it corresponds to the estimated relative frequency under the assumption of independence between Y and Z.
} 


\value{

The required graphical representation is drawn using standard \pkg{graphics} facilities.

}

\author{
 Marcello D'Orazio \email{mdo.statmatch@gmail.com} 
}

\seealso{
\code{\link[StatMatch]{Frechet.bounds.cat}}
}

\examples{

# # compute the tables required by Frechet.bounds.cat()
# freq.xA <- xtabs(~sex+c.age, data=samp.A)
# freq.xB <- xtabs(~sex+c.age, data=samp.B)
# freq.xy <- xtabs(~sex+c.age+c.neti, data=samp.A)
# freq.xz <- xtabs(~sex+c.age+labour5, data=samp.B)
# 
# # apply Frechet.bounds.cat()
# bounds.yz <- Frechet.bounds.cat(tab.x=freq.xA+freq.xB, tab.xy=freq.xy,
#                                 tab.xz=freq.xz, print.f="data.frame")
# 
# 
# plot.bounds(bounds.yz)


}

\keyword{multivariate}
\keyword{dplot}