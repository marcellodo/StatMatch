\name{comp.cont}
\alias{comp.cont}
\title{Empirical comparison of two estimated distributions of the same continuous variable}

\description{
  This function estimates the \dQuote{closeness} of distributions of the same continuous variable(s) but estimated from different data sources.
}

\usage{
comp.cont(data.A, data.B, xlab.A, xlab.B = NULL, w.A = NULL, w.B = NULL, ref = FALSE) 
}

\arguments{

\item{data.A}{
 A dataframe or matrix containing the variable of interest \code{xlab.A} and eventual survey weights \code{w.A}. 
}

\item{data.B}{ 
	A dataframe or matrix containing the variable of interest \code{xlab.B} and eventual associated survey weights \code{w.B}. 
}

\item{xlab.A}{
 Character string providing the name of the variable in \code{data.A} whose estimated distribution should be compared with that estimated from \code{data.B}.
}

\item{xlab.B}{
 Character string providing the name of the variable in \code{data.B} whose distribution should be  compared with that estimated from \code{data.A}. If \code{xlab.B=NULL} (default) then it assumed \code{xlab.B=xlab.A}.
}

\item{w.A}{
   Character string providing the name of the optional weighting variable in \code{data.A} that, in case, should be used to estimate the distribution of \code{xlab.A}
}

\item{w.B}{
   Character string providing the name of the optional weighting variable in \code{data.B} that, in case, should be used to estimate the distribution of \code{xlab.B}
}

\item{ref}{
	Logical.  When \code{ref = TRUE}, the distribution of \code{xlab.B} estimated from \code{data.B} is considered the reference distribution (true or reliable estimate of distribution). Affects some estimation procedures as explained in the Details. 
}

}


\details{
As a first output, the function returns some well--known summary measures (min, Q1, median, mean, Q3, max and sd) estimated from the available input data sources. 

Secondly this function performs a comparison between the  quantiles estimated from \code{data.A} and \code{data.B}; in particular, the average of the absolute value of the differences as well as the average of the squared differences are returned. The number of estimated percentiles depends on the minimum between the two sample sizes. Only quartiles are calculated when min(n.A, n.B)<20, deciles are estimated when min(n.A, n.B)>=20 and min(n.A, n.B)<=30, finally quantiles for \code{probs=seq(from = 0.05,to = 0.95,by = 0.05)} are estimated when min(n.A, n.B)>30. When the survey weights are available (indicated with th arguments \code{w.A} and/or \code{w.B}) they are used in estimating the quantiles by calling the function \code{\link[Hmisc]{wtd.quantile}} in the package \pkg{Hmisc}.
  
The function estimates also the dissimilarities between the estimated empirical distribution function. The measures considered are the maximum of the absolute differences, the sum between the maximum differences inverting the terms in the difference and the average of the absolute value of the differences. When the weights are provided they are used in estimating the empirical cumulative distribution function. Note that when \code{ref=TRUE} the estimation of the density and of the empirical cumulative distribution are guided by the data in \code{data.B}.

The final output is the total variation distance, the overlap and the Hellinger distance calculated considering the transformed categorized variable. The breaks to categorize the variable are decided according to the Freedman-Diaconis rule (\code{\link[grDevices]{nclass}}) and, in this case, when \code{ref=TRUE} the IQR is estimated solely on \code{data.B}, whereas with \code{ref=FALSE} it is estimated by joining the two data sources. 
When present, the weights are used in estimating the relative frequencies of the categorized variable. 
For additional details on these distances please see (\code{\link[StatMatch]{comp.prop}})

}

\value{
A \code{list} object with four components. 

\item{summary}{A matrix with summaries of \code{xlab.A} estimated on \code{data.A} and summaries of \code{xlab.B} estimated on \code{data.B}}

\item{diff.Qs}{ Average of absolute and squared differences between the quantiles of \code{xlab.A} estimated on \code{data.A} and the corresponding ones of \code{xlab.B} estimated on \code{data.B}}

\item{dist.ecdf}{Dissimilarity measures between the estimated empirical cumulative distribution functions.}

\item{dist.discr}{Distance between the distributions after discretization of the target variable.}

}

\references{

Bellhouse D.R. and J. E. Stafford (1999). \dQuote{Density Estimation from Complex Surveys}. \emph{Statistica Sinica}, \bold{9}, 407--424.
}

\author{
 Marcello D'Orazio \email{mdo.statmatch@gmail.com} 
}

\seealso{ 
\code{\link[StatMatch]{plotCont}}, \code{\link[StatMatch]{comp.prop}} 
}

\examples{
data(samp.A)
data(samp.B)

comp.cont(data.A = samp.A, data.B = samp.B, xlab.A = "age")

comp.cont(data.A = samp.A, data.B = samp.B, xlab.A = "age",
          w.A = "ww", w.B = "ww")

}

\keyword{multivariate}
