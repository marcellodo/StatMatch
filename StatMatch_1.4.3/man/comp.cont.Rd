\name{comp.cont}
\alias{comp.cont}
\title{Compares two distributions of the same continuous variable}

\description{
  This function estimate the \dQuote{closeness} of the distributions of the same continuous variable(s) but estimated from different data sources.
}

\usage{
comp.cont(data.A, data.B, xlab.A, xlab.B = NULL, w.A = NULL, 
          w.B = NULL, ref = FALSE) 
}

\arguments{

\item{data.A}{
 A dataframe or matrix containing the variable of interest \code{xlab.A} and eventual associated survey weights \code{w.A}. 
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
 Logical.  When \code{ref = TRUE}, the distribution of \code{xlab.B} estimated from \code{data.B} is considered the reference distribution (true or reliable estimate of distribution).
}

}


\details{
This function computes well--known summary measures (min, Q1, median, mean, Q3, max and sd) estimated from the available data. In addition, quantiles estimated \code{data.A} are compared with those estimated from \code{data.B}, and the average of the absolute value of the differences as well as the average of the squared differences are returned. Note that the number of estimated percentiles depends on the minimum between the two sample sizes. Note that the number of estimated percentiles depends on the minimum between the two sample sizes. Only quartiles are calculated when min(n.A, n.B)<=50; quintiles are estimated when min(n.A, n.B)>50 and min(n.A, n.B)<=150; deciles are estimated when min(n.A, n.B)>150 and min(n.A, n.B)<=250; finally quantiles for \code{probs=seq(from = 0.05,to = 0.95,by = 0.05)} are estimated when min(n.A, n.B)>250. When survey weights are available (indicated through \code{w.A} and/or \code{w.B}) they are used in estimating the quantiles by calling the function \code{\link[Hmisc]{wtd.quantile}} in the package \pkg{Hmisc}.
  
Dissimilarities between the estimated empirical distribution function are calculated. The measures considered are the maxim value of the differences, the sum between the absolute value of both minimum and maximum and the average of the absolute differences. When weights are provided they are used in estimating the empirical cumulative distribution function. Note that when \code{ref=TRUE} the estimation of the density and the empirical cumulative distribution are guided by the data in  \code{data.B}.

Finally is calculated the total variation distance, the overlap and the Hellinger, on the transformed categorized variable. Note that the breaks to categorize the variable are decided according to the Freedman-Diaconis rule (\code{\link[grDevices]{nclass}}) and, in this case, with \code{ref=TRUE} the IQR is estimated solely on \code{data.A}, whereas with \code{ref=FALSE} it estimated by joining the two data sources. 
When present, the weights are used in estimating the relative frequencies of the categorized variable. 

\emph{total variation distance}: 

\deqn{\Delta_{AB} = \frac{1}{2} \sum_{j=1}^J \left| p_{A,j} - p_{B,j}  \right|}{D = (1/2) * sum_j |p_A,j - p_B,j| }

where \eqn{p_{s,j}}{p_s,j} are the relative frequencies (\eqn{0 \leq p_{s,j} \leq 1}{0 <= p_s,j <= 1}). The dissimilarity index ranges from 0 (minimum dissimilarity) to 1. The total variation distance comes along with its complement to 1, said \dQuote{overlap} between distributions.
 

the \emph{Hellinger's distance}:

\deqn{d_{H,AB} = \sqrt{  \frac{1}{2} \sum_{j=1}^J \left( \sqrt{p_{A,j}} - \sqrt{p_{B,j}}  \right)^2 } }{d_H = ((1/2) * sum_j (sqrt(p_A,j) - sqrt(p_B,j))^(1/2))^(1/2) }

It is a dissimilarity measure ranging from 0 (distributions are equal) to 1 (max dissimilarity).  It satisfies all the properties of a distance measure (\eqn{0 \leq d_{H,AB} \leq 1}{0 <= d_H <= 1};  symmetry and triangle inequality).  
Hellinger's distance is related to the total variation distance, and it is possible to show that:

\deqn{d_{H,AB}^2 \leq  \Delta_{AB} \leq d_{H,AB}\sqrt{2} }{ d_H^2 <= D <= d_H * sqrt(2) }

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
\code{\link[StatMatch]{plotCont}}
}

\keyword{multivariate}
