% --- Source file: fact2dummy.Rd ---
\name{plotCont}
\alias{plotCont}
\title{graphical comparison of the estimated distributions for the same continuous variable.}

\description{
Compares graphically the estimated distributions for the same continuous variable using data coming from two different data sources.
}

\usage{
plotCont(data.A, data.B, xlab.A, xlab.B=NULL, w.A=NULL, w.B=NULL,
         type="density", ref=FALSE) 
}

\arguments{

\item{data.A}{
  	A dataframe or matrix containing the variable of interest \code{xlab.A} and eventual associated survey weights \code{w.A}. 
}

\item{data.B}{
  	A dataframe or matrix containing the variable of interest \code{xlab.B} and eventual associated survey weights \code{w.B}. 
}

\item{xlab.A}{
  	Character string providing the name of the variable in \code{data.A} whose distribution should be represented graphically and compared with that estimated from \code{data.B}.
}

\item{xlab.B}{
  	Character string providing the name of the variable in \code{data.B} whose distribution should be represented graphically and compared with that estimated from \code{data.A}. If \code{xlab.B=NULL} (default) then it assumed \code{xlab.B=xlab.A}.
}

\item{w.A}{
   Character string providing the name of the optional weighting variable in \code{data.A} that, in case, should be used to estimate the distribution of \code{xlab.A}
}

\item{w.B}{
   Character string providing the name of the optional weighting variable in \code{data.B} that, in case, should be used to estimate the distribution of \code{xlab.B}
}

\item{type}{
   A character string indicating the type of graphical representation that should be used to compare the estimated distributions of \code{xlab.A} and \code{xlab.B}. By default (\code{type="density"}) density plots are used. Other possible options are \dQuote{ecdf}, \dQuote{qqplot}, \dQuote{qqshift} and \dQuote{hist}. See Details for more information.
}

\item{ref}{
   Logical, indicating whether the distribution estimated from \code{data.B} should be considered the reference or not. Default \code{ref=FALSE}. when Default \code{ref=TRUE} the estimation of the histograms, the density and the empirical cumulative distribution function are guided by data in  \code{data.B} 
}

}
  
\details{
This function compares graphically the distribution of the same variable but estimated from data coming from two different data sources. The graphical comparison con be done in different manners. 
When \code{type="hist"} the continuous variable is categorized and the corresponding histograms, estimated from \code{data.A} and \code{data.B}, are compared. When present, the weights are used in estimating the relative frequencies. Note that the breaks to categorize the variable are decided according to the Freedman-Diaconis rule (\code{\link[grDevices]{nclass}}) and, in this case, with \code{ref=TRUE} the IQR is estimated solely on \code{data.B}, whereas with \code{ref=FALSE} it is estimated by joining the two data sources.
 
With \code{type="density"} the density plots are drawn; when available the weights are used in the estimation of the density that are based on the histograms (as suggested by Bellhouse and Stafford, 1999). When\code{type="ecdf"} the comparison relies on the empirical cumulative distribution function, that can be estimated considering the weights. Note that when \code{ref=TRUE} the estimation of the density and the empirical cumulative distribution are guided by the data in  \code{data.B}.

The comparison is based on percentiles with \code{type="qqplot"} and \code{type="qqshift"}. In the first case, the function draws a scatterplot (red dots) of the estimated percentiles of \code{xlab.A} vs. those of \code{xlab.B}; the dashed line indicated the ideal situation of equality of percentiles (points lying on the line). When \code{type="qqshift"} the scatterplot refers to (percentiles.A - percentiles.B) vs. percentiles.B; in this case the points lying on horizontal line passing through 0 indicate equality (difference equal to 0). Note that the number of estimated percentiles depends on the minimum between the two sample sizes. Only quartiles are calculated when min(n.A, n.B)<=50; quintiles are estimated when min(n.A, n.B)>50 and min(n.A, n.B)<=150; deciles are estimated when min(n.A, n.B)>150 and min(n.A, n.B)<=250; finally quantiles for \code{probs=seq(from = 0.05,to = 0.95,by = 0.05)} are estimated when min(n.A, n.B)>250. When survey weights are available (indicated through \code{w.A} and/or \code{w.B}) they are used in estimating the quantiles by calling the function \code{\link[Hmisc]{wtd.quantile}} in the package \pkg{Hmisc}.
} 


\value{

The required graphical representation is drawn using the \pkg{ggplot2} facilities.

}

\references{

Bellhouse D.R. and J. E. Stafford (1999). \dQuote{Density Estimation from Complex Surveys}. \emph{Statistica Sinica}, \bold{9}, 407--424.

}


\author{
 Marcello D'Orazio \email{mdo.statmatch@gmail.com} 
}

\seealso{
\code{\link[StatMatch]{comp.cont}}
}

\examples{

# plotCont(data.A = samp.A, data.B = samp.B, xlab.A="age")
# plotCont(data.A = samp.A, data.B = samp.B, xlab.A="age", w.A = "ww")

}

\keyword{multivariate}
\keyword{dplot}