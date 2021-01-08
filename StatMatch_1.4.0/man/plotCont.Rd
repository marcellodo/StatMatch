% --- Source file: fact2dummy.Rd ---
\name{plotCont}
\alias{plotCont}
\title{graphical comparison of the estimated distributions for the same continuous variable.}

\description{
Compares graphically the estimated distributions for the same continuous variable using data coming from two different data sources.
}

\usage{
plotCont(data.A, data.B, xlab.A, xlab.B=NULL, w.A=NULL, w.B=NULL,
         type="density") 
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
   A character string indicating the type of graphical output that should be used to compare the estimated distributions of \code{xlab.A} and \code{xlab.B}. By default (\code{type="density"}) density plots are used. Other possible options are \dQuote{qqplot}, \dQuote{qqshift} and \dQuote{hist}. See Details for more information.
}

}
  
\details{
This function compares graphically distribution of the same variable but estimated from data coming from two different data sources. The graphical comparison con be done in four different manners. When \code{type="density"} the density plots are drawn; when available the weights are used in the estimation of the density. 

The comparison is based on percentiles with \code{type="qqplot"} and \code{type="qqshift"}. In the first case, the function draws a scatterplot (red dots) of the estimated percentiles of \code{xlab.A} vs. those of \code{xlab.B}; the dashed line indicated the ideal situation of equality of percentiles (points lying on the line). When \code{type="qqshift"} the scatterplot refers to (percentiles.A - percentiles.B) vs. percentiles.A; in this case the points lying on horizontal line passing through 0 indicate equality (difference equal to 0). Note that the number of estimated percentiles depends on the minimum between the two sample sizes. Only quartiles are calculated when min(n.A, n.B)<20, deciles are estimated when min(n.A, n.B)>=20 and min(n.A, n.B)<=30, finally quantiles for \code{probs=seq(from = 0.05,to = 0.95,by = 0.05)} are estimated when min(n.A, n.B)>30. When survey weights are available (indicated through \code{w.A} and/or \code{w.B}) they are used in estimating the quantiles by calling the function \code{\link[Hmisc]{wtd.quantile}} in the package \pkg{Hmisc}.

Finally, when \code{type="hist"} the continuous variable is categorized and the corresponding histograms, estimated from \code{data.A} and \code{data.B}, are compared. Also in this case, when present, the weights are used in estimating the relative frequencies.  
   
} 


\value{

The required graphical representation is drawn using the \pkg{ggplot2} facilities.

}

\author{
 Marcello D'Orazio \email{mdo.statmatch@gmail.com} 
}


\examples{

# plotCont(data.A = samp.A, data.B = samp.B, xlab.A="age")
# plotCont(data.A = samp.A, data.B = samp.B, xlab.A="age", w.A = "ww")

}

\keyword{multivariate}
\keyword{dplot}