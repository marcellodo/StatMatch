% --- Source file:  ---
\name{Fbwidths.by.x}
\alias{Fbwidths.by.x}
\title{Computes the Frechet bounds of cells in a contingency table by considering all the possible subsets of the common variables.}

\description{
This function permits to compute the bounds for cell probabilities in the contingency table Y vs. Z starting from the marginal tables (\bold{X} vs. Y), (\bold{X} vs. Z) and the joint distribution of the \bold{X} variables, by considering all the possible subsets of the \bold{X} variables.  In this manner it is possible to identify which subset of the \bold{X} variables produces the major reduction of the average width of conditional bounds.
}

\usage{
Fbwidths.by.x(tab.x, tab.xy, tab.xz, deal.sparse="discard", 
          nA=NULL, nB=NULL, ...)
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

%\item{compress.sum}{
%Logical (default \code{FALSE)}. If \code{TRUE} reduces the information saved in \code{sum.unc}. See Value for further information.
%}
\item{deal.sparse}{
Text, how to estimate the cell relative frequencies when dealing with too sparse tables. When \code{deal.sparse="discard"} (default) no estimation is performed if \code{tab.xy} or \code{tab.xz} is too sparse. When \code{deal.sparse="relfreq"} the standard estimator (cell count divided by the sample size) is considered. 
Note that here sparseness is measured by number of cells with respect to the sample size; sparse table are those where the number of cells exceeds  the sample size (see Details).
}

\item{nA}{
Integer, sample size of file A used to estimate \code{tab.xy}. If  \code{NULL}, it is obtained as sum of frequencies in\code{tab.xy}.
}

\item{nB}{
Integer, sample size of file B used to estimate \code{tab.xz}. If  \code{NULL}, it is obtained as sum of frequencies in\code{tab.xz}.
}

\item{...}{
Additional arguments that may be required when deriving an estimate of uncertainty by calling \code{\link[StatMatch]{Frechet.bounds.cat}}.
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

Provided that uncertainty, measured in terms of \eqn{\bar{d}}{av(d)}, tends to reduce when conditioning on a higher number of \bold{X} variables. Two penalties are introduced to account for the additional number of cells to be estimated when adding a X variable. The first penalty, introduced in D'Orazio et al. (2017), is:

\deqn{g_1=log\left( 1 + \frac{H_{D_m}}{H_{D_Q}} \right) }{g1 = log(1 + H_Dm/H_DQ )}
    
Where \eqn{H_{D_m}}{H_Dm} is the number of cell in the table obtained by crossing the given subset of \bold{X} variables and the \eqn{H_{D_Q}}{H_DQ} is the number of cell in the table achieved by crossing all the available \bold{X} variables. 
A second penalty takes into account the number of cells to estimate with respect to the sample size (D'Orazio et al., 2019). It is obtained as:

\deqn{g_2 = max \left[ \frac{1}{n_A - H_{D_m} \times J}, \frac{1}{n_B - H_{D_m} \times K}  \right]}{ g2=max[1/(nA - H_Dm*J), 1/(nB - H_Dm*K)]}

with \eqn{n_A > H_{D_m} \times J}{nA > H_Dm*J}  and  \eqn{n_B > H_{D_m} \times K}{nB > H_Dm*K}. In practice, it is considered the number of cells to estimate compared to the sample size. This criterion is considered to measure sparseness too. In particular, for the purposes of this function, tables are NOT considered sparse when: 

\deqn{min\left[ \frac{n_A}{H_{D_m} \times J}, \frac{n_B}{H_{D_m} \times K} \right] > 1 }{ min[ nA/(H_Dm*J), nB/(H_Dm*K) ] > 1}

This rule is applied when deciding how to proceed with estimation in case of sparse table (argument \code{deal.sparse}). 
Note that sparseness can be measured in different manners. The outputs include also the empty cells in each table (due to statistical zeros or structural zeros) and the Cohen's effect size with respect to the case of uniform distribution of frequencies across cells (the value 1/no.of.cells in every cell):

\deqn{\omega_{eq} = \sqrt{H \sum_{h=1}^{H} (\hat{p}_h - 1/H)^2 } }{ w_eq = ( H*sum((p_h - 1/H)^2))^(1/2) }

values of \eqn{\omega_{eq}}{w_eq>2} jointly with  \eqn{n/H \leq 1}{n/H<=1} usually indicate severe sparseness. 
}

\value{

A list with the estimated bounds for the cells in the table of Y vs. Z for each possible subset of the \bold{X} variables.  The final component in the list, \code{sum.unc}, is a data.frame that summarizes the main results. In particular, it reports the number of \bold{X} variables (\code{"x.vars"}), the number of cells in each of the input tables and the cells with frequency equal to 0 (columns ending with \code{freq0} ). Moreover, it reported the value (\code{"av.n"}) of the rule used to decide whether we are dealing with a sparse case (see Details) and the Cohen's effect size measured for the table crossing the considered combination of the X variables. 
Finally, it is provided the average width of the uncertainty intervals (\code{"av.width"}), the penalty terms g1 and g2 (\code{"penalty1"} and \code{"penalty2"} respectively), and the penalized average widths (\code{"av.width.pen1"} and \code{"av.width.pen2"}, where av.width.pen1=av.width+pen1 and av.width.pen2=av.width+pen2).

%When \code{compress.sum = TRUE} the data.frame \code{sum.unc} will show a combination of the \bold{X} variables only if it determines a reduction of the (\code{"av.width"}) when compared to the preceding one.

%Note that in the presence of too many cells with 0s in the input contingency tables is an indication of sparseness; this is an unappealing situation when estimating the cells' relative frequencies needed to derive the bounds; in such cases the corresponding results may be unreliable. A possible alternative way of working consists in estimating the required parameters by considering a pseudo-Bayes estimator (see \code{\link[StatMatch]{pBayes}}); in practice the input \code{tab.x}, \code{tab.xy} and \code{tab.xz} should be the ones provided by the \code{\link[StatMatch]{pBayes}} function.

}

\references{

%Conti P.L, Marella, D., Scanu, M. (2012) \dQuote{Uncertainty Analysis in Statistical Matching}. \emph{Journal of Official Statistics}, \bold{28}, pp. 69--88.

D'Orazio, M., Di Zio, M. and Scanu, M. (2006). \emph{Statistical Matching: Theory and Practice.} Wiley, Chichester.

D'Orazio, M., Di Zio, M. and Scanu, M. (2017). ``The use of uncertainty to choose matching variables in statistical matching''. \emph{International Journal of Approximate Reasoning }, 90, pp. 433-440.

D'Orazio, M., Di Zio, M. and Scanu, M. (2019). ``Auxiliary variable selection in a a statistical matching problem''. In Zhang, L.-C. and Chambers, R. L. (eds.) \emph{Analysis of Integrated Data}, Chapman & Hall/CRC (Forthcoming).
}


\author{
 Marcello D'Orazio \email{mdo.statmatch@gmail.com} 
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
suppressWarnings(RNGversion("3.5.0"))
set.seed(4567)
lab.A <- sample(nrow(quine), 70, replace=TRUE)
quine.A <- quine[lab.A, 1:4]
quine.B <- quine[-lab.A, c(1:3,6)]

# compute the tables required by Fbwidths.by.x()
freq.xA <- xtabs(~Eth+Sex+Age, data=quine.A)
freq.xB <- xtabs(~Eth+Sex+Age, data=quine.B)

freq.xy <- xtabs(~Eth+Sex+Age+Lrn, data=quine.A)
freq.xz <- xtabs(~Eth+Sex+Age+c.Days, data=quine.B)

# apply Fbwidths.by.x()
bounds.yz <- Fbwidths.by.x(tab.x=freq.xA+freq.xB, tab.xy=freq.xy,
        tab.xz=freq.xz)

bounds.yz$sum.unc


}

\keyword{multivariate}