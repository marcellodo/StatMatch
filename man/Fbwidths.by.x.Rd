% --- Source file:  ---
\name{Fbwidths.by.x}
\alias{Fbwidths.by.x}
\title{Computes the Frechet bounds of cells in a contingency table by considering all the possible subsets of the common variables.}

\description{
This function permits to compute the bounds for cell probabilities in the contingency table Y vs. Z starting from the marginal tables (\bold{X} vs. Y), (\bold{X} vs. Z) and the joint distribution of the \bold{X} variables, by considering all the possible subsets of the \bold{X} variables.  In this manner it is possible to identify which subset of the \bold{X} variables produces the major reduction of the uncertainty.
}

\usage{
Fbwidths.by.x(tab.x, tab.xy, tab.xz, compress.sum=FALSE)
}

\arguments{

\item{tab.x}{
A \R table crossing the \bold{X} variables.  This table must be obtained by using the function \code{\link[stats]{xtabs}} or \code{\link[base]{table}}, e.g. \cr
\code{tab.x <- xtabs(~x1+x2+x3, data=data.all)}.
}

\item{tab.xy}{
A \R table of \bold{X} vs. Y variable.  This table must be obtained by using the function \code{\link[stats]{xtabs}} or \code{\link[base]{table}}, e.g. \cr
\code{table.xy <- xtabs(~x1+x2+x3+y, data=data.A)}.

A single categorical Y variables is allowed.  One or more categorical variables can be considered as \bold{X} variables (common variables).  The same \bold{X} variables in \code{tab.x} must be available in \code{tab.xy}.  Moreover, it is assumed that the joint distribution of the \bold{X} variables computed from \code{tab.xy} is equal to \code{tab.x}; a warning is produced if this is not true.
}

\item{tab.xz}{
A \R table of \bold{X} vs. Z variable.  This table must be obtained by using the function \code{\link[stats]{xtabs}} or \code{\link[base]{table}}, e.g. \cr
\code{tab.xz <- xtabs(~x1+x2+x3+z, data=data.B)}.

A single categorical Z variable is allowed.  One or more categorical variables can be considered as \bold{X} variables (common variables).  The same \bold{X} variables in \code{tab.x} must be available in \code{tab.xz}.  Moreover, it is assumed that the joint distribution of the \bold{X} variables computed from \code{tab.xz} is equal to \code{tab.x}; a warning is produced if this is not true.
}

\item{compress.sum}{
Logical (default \code{FALSE)}. If \code{TRUE} reduces the information saved in \code{sum.unc}. See Value for further information.
}


}

\details{
This function permits to compute the Frechet bounds for the frequencies in the contingency table of Y vs. Z, starting from the conditional distributions P(Y|\bold{X}) and P(Z|\bold{X}) (for details see \cr
\code{\link[StatMatch]{Frechet.bounds.cat}}), by considering all the possible subsets of the \bold{X} variables.  In this manner it is possible to identify the subset of the \bold{X} variables, with highest association with both Y and Z, that permits to reduce the uncertainty concerning the distribution of Y vs. Z. 

%The overall uncertainty is measured by considering the suggestion in Conti \emph{et al.} (2012):

%\deqn{ \hat{\Delta} = \sum_{i,j,k} ( p^{(up)}_{Y=j,Z=k} - p^{(low)}_{Y=j,Z=k} ) \times p_{Y=j|X=i} \times p_{Z=k|X=i} \times p_{X=i}  
%}{ D = sum_(i,j,k) ( p^(up)_(Y=j,Z=k) - p^(low)_(Y=j,Z=k) ) * p_(Y=j|X=i) * p_(Z=k|X=i) * p_(X=i) }


The uncertainty is measured by the average of the widths of the bounds for the cells in the table Y vs. Z:

\deqn{ \bar{d} = \frac{1}{J \times K} \sum_{j,k} ( p^{(up)}_{Y=j,Z=k} - p^{(low)}_{Y=j,Z=k} )}{d=(1/(J*K))*sum_(j,k)(p^(up)_(Y=j,Z=k) - p^(low)_(Y=j,Z=k))}

For details see \code{\link[StatMatch]{Frechet.bounds.cat}}.

}


\value{

A list with the estimated bounds for the cells in the table of Y vs. Z for each possible subset of the \bold{X} variables.  The final component  \code{sum.unc} is a data.frame that summarizes the main findings. In particular it reports the number of \bold{X} variables (\code{"x.vars"}), the number of cells in the each of the input tables and the corresponding number of cells with frequency equal to 0 (columns ending with \code{freq0} ). Then it is provided the average width of the uncertainty intervals (\code{"av.width"}) and its relative value (\code{"rel.av.width"}) when compared with the average widths of the uncertainty intervals when no \bold{X} variables are considered (i.e. \code{unconditioned} \code{"av.width"}, reported in the first row of the data.frame).

When \code{compress.sum = TRUE} the data.frame \code{sum.unc} will show a combination of the \bold{X} variables only if it determines a reduction of the (\code{"av.width"}) when compared to the preceding one.

Note that in the presence of too many cells with 0s in the input contingency tables is an indication of sparseness; this is an unappealing situation when estimating the cells' relative frequencies needed to derive the bounds; in such cases the corresponding results may be unreliable. A possible alternative way of working consists in estimating the required parameters by considering a pseudo-Bayes estimator (see \code{\link[StatMatch]{pBayes}}); in practice the input \code{tab.x}, \code{tab.xy} and \code{tab.xz} should be the ones provided by the \code{\link[StatMatch]{pBayes}} function.



}

\references{

Ballin, M., D'Orazio, M., Di Zio, M., Scanu, M. and Torelli, N. (2009) \dQuote{Statistical Matching of Two Surveys with a Common Subset}. \emph{Working Paper}, \bold{124}. Dip. Scienze Economiche e Statistiche, Univ. di Trieste, Trieste. 

%Conti P.L, Marella, D., Scanu, M. (2012) \dQuote{Uncertainty Analysis in Statistical Matching}. \emph{Journal of Official Statistics}, \bold{28}, pp. 69--88.

D'Orazio, M., Di Zio, M. and Scanu, M. (2006). \emph{Statistical Matching: Theory and Practice.} Wiley, Chichester.

}


\author{
 Marcello D'Orazio \email{madorazi@istat.it} 
}

\seealso{ 
\code{\link[StatMatch]{Frechet.bounds.cat}}, \code{\link[StatMatch]{harmonize.x}}
}

\examples{

data(quine, package="MASS") #loads quine from MASS
str(quine)
quine$c.Days <- cut(quine$Days, c(-1, seq(0,50,10),100))
table(quine$c.Days)


# split quine in two subsets
set.seed(4567)
lab.A <- sample(nrow(quine), 70, replace=TRUE)
quine.A <- quine[lab.A, 1:4]
quine.B <- quine[-lab.A, c(1:3,6)]

# compute the tables required by Fbwidths.by.x()
freq.x <- xtabs(~Eth+Sex+Age, data=quine.A)
freq.xy <- xtabs(~Eth+Sex+Age+Lrn, data=quine.A)
freq.xz <- xtabs(~Eth+Sex+Age+c.Days, data=quine.B)

# apply Fbwidths.by.x()
bounds.yz <- Fbwidths.by.x(tab.x=freq.x, tab.xy=freq.xy,
        tab.xz=freq.xz)

bounds.yz$sum.unc

# input tables estimated with pBayes()

pf.x <- pBayes(x=freq.x)
pf.xy <- pBayes(x=freq.xy)
pf.xz <- pBayes(x=freq.xz)

bounds.yz.p <- Fbwidths.by.x(tab.x = pf.x$pseudoB, 
							 tab.xy = pf.xy$pseudoB,
							 tab.xz = pf.xz$pseudoB)

}

\keyword{multivariate}