% --- Source file:  ---
\name{selMtc.by.unc}
\alias{selMtc.by.unc}
\title{Identifies the best combination if matching variables in reducing uncertainty in estimation the contingency table Y vs. Z.}

\description{
This function identifies the ``best'' subset of matching variables in terms of reduction of uncertainty when estimating relative frequencies in the contingency table Y vs. Z. The sequential procedure presented in D'Orazio \emph{et al.}  (2017 and 2019) is implemented. This procedure avoids exploring all the possible combinations of the available X variables as in \code{\link[StatMatch]{Fbwidths.by.x}}.
}

\usage{
selMtc.by.unc(tab.x, tab.xy, tab.xz, corr.d=2, 
                    nA=NULL, nB=NULL, align.margins=FALSE) 
}

\arguments{

\item{tab.x}{
A \R table crossing the \bold{X} variables.  This table must be obtained by using the function \code{\link[stats]{xtabs}} or \code{\link[base]{table}}, e.g. \cr
\code{tab.x <- xtabs(~x1+x2+x3, data=data.all)}. \strong{A minimum number of 3 variables is needed}.
}

\item{tab.xy}{
A \R table of \bold{X} vs. Y variable.  This table must be obtained by using the function \code{\link[stats]{xtabs}} or \code{\link[base]{table}}, e.g. \cr
\code{table.xy <- xtabs(~x1+x2+x3+y, data=data.A)}.

A single categorical Y variables is allowed.  At least \strong{three} categorical variables should be considered as \bold{X} variables (common variables).  The same \bold{X} variables in \code{tab.x} must be available in \code{tab.xy}.   Usually, it is assumed that the joint distribution of the \bold{X} variables computed from \code{tab.xy} is equal to \code{tab.x} (a warning appears if any absolute difference is greater than \code{tol}). Note that when the marginal distribution of  \bold{X} in \code{tab.xy} is not equal to that of \code{tab.x} it is possible to align them before computations (see argument \code{align.margins}).

}

\item{tab.xz}{
A \R table of \bold{X} vs. Z variable.  This table must be obtained by using the function \code{\link[stats]{xtabs}} or \code{\link[base]{table}}, e.g. \cr
\code{tab.xz <- xtabs(~x1+x2+x3+z, data=data.B)}.

A single categorical Z variable is allowed.  At least \strong{three} categorical variables should be considered as \bold{X} variables (common variables).  The same \bold{X} variables in \code{tab.x} must be available in \code{tab.xz}.  Usually, it is assumed that the joint distribution of the \bold{X} variables computed from \code{tab.xz} is equal to \code{tab.x} (a warning appears if any absolute difference is greater than \code{tol}). Note that when the marginal distribution of  \bold{X} in \code{tab.xz} is not equal to that of \code{tab.x} it is possible to align them before computations (see argument \code{align.margins}).
}

\item{corr.d}{
Integer, indicates the penalty that should be introduced in estimating the uncertainty by means of the average width of cell bounds. When \code{corr.d=1}  the penalty being considered is the one introduced in D'Orazio \emph{et al.} (2017) (i.e. penalty1 in \code{\link[StatMatch]{Fbwidths.by.x}}). When \code{corr.d=2} (default) it is considered a penalty suggested in D'Orazio \emph{et al.} (2019) (indicated as ``penalty2''  in   \code{\link[StatMatch]{Fbwidths.by.x}}). Finally, no penalties are considered when \code{corr.d=0}.
}

\item{nA}{
Integer, sample size of file A used to estimate \code{tab.xy}. If  \code{NULL} is obtained as sum of frequencies in \code{tab.xy}.
}

\item{nB}{
Integer, sample size of file B used to estimate \code{tab.xz}. If  \code{NULL} is obtained as sum of frequencies in \code{tab.xz}.
}

\item{align.margins}{
Logical (default \code{FALSE}). When when \code{TRUE} the distribution of \bold{X} variables in \code{tab.xy} is aligned with the distribution resulting from \code{tab.x}, without affecting the marginal distribution of Y. Similarly the distribution of \bold{X} variables in \code{tab.xz} is aligned with the distribution resulting from \code{tab.x}, without affecting the marginal distribution of Z. The alignment is performed by running IPF algorithm as implemented in the function \code{\link[mipfp]{Estimate}} in the package \pkg{mipfp}. To avoid lack of convergence due to combinations of Xs encountered in one table but not in the other (statistical 0s), before running IPF a small constant (1e-06) is added to empty cells in \code{tab.xy} and \code{tab.xz}.
}

}

\details{
This function follows the sequential procedure described in D'Orazio \emph{et al.} (2017, 2019) to identify the combination of common variables most effective in reducing uncertainty when estimating the contingency table Y vs. Z. Initially, the available Xs are ordered according to the reduction of average width of uncertainty bounds when conditioning on each of them. Then in each step one the remaining X variables is added until the table became too sparse; in practice the procedure stops when: 

\deqn{min\left[ \frac{n_A}{H_{D_m} \times J}, \frac{n_B}{H_{D_m} \times K} \right] \leq 1 }{ min[ nA/(H_Dm*J), nB/(H_Dm*K) ] <= 1} 

For major details see also \code{\link[StatMatch]{Fbwidths.by.x}}.

}

\value{

A list with the main outcomes of the procedure.

\item{ini.ord}{
Average width of uncertainty bounds when conditioning on each of the available X variables. Variable most effective in reducing uncertainty comes first. The ordering determines the order in which they are entered in the sequential procedure. 
}

\item{list.xs}{
List with the various combinations of the matching variables being considered in each step.
}

\item{av.df}{
Data.frame with all the relevant information for each of combination of X variables. The last row corresponds to the combination of the X variables identified as the best in reducing average width of uncertainty bounds (penalized or not depending on the input argument \code{corr.d}). For each combination of X variables the following additional information are reported:  the number of cells (name starts with ``\code{nc}'');  the number of empty cells (name starts with ``\code{nc0}''; the average relative frequency (name starts with ``\code{av.crf}''); sparseness measured as Cohen's effect size with respect to equiprobability (uniform distribution across cells). Finally there are the value of the stopping criterion (``\code{min.av}''), the unconditioned average width of uncertainty bounds (``\code{avw}''), the penalty term  (``\code{penalty}'') and the penalized width  (``\code{avw.pen}''; \code{avw.pen=avw+penalty}).
}

}

\references{

D'Orazio, M., Di Zio, M. and Scanu, M. (2006). \emph{Statistical Matching: Theory and Practice.} Wiley, Chichester.

D'Orazio, M., Di Zio, M. and Scanu, M. (2017). ``The use of uncertainty to choose matching variables in statistical matching''. \emph{International Journal of Approximate Reasoning}, 90, pp. 433-440.

D'Orazio, M., Di Zio, M. and Scanu, M. (2019). ``Auxiliary variable selection in a a statistical matching problem''. In Zhang, L.-C. and Chambers, R. L. (eds.) \emph{Analysis of Integrated Data}, Chapman & Hall/CRC (forthcoming).
}


\author{
 Marcello D'Orazio \email{mdo.statmatch@gmail.com} 
}

\seealso{ 
\code{\link[StatMatch]{Fbwidths.by.x}}, \code{\link[StatMatch]{Frechet.bounds.cat}}
}

\examples{

data(quine, package="MASS") #loads quine from MASS
str(quine)
quine$c.Days <- cut(quine$Days, c(-1, seq(0,50,10),100))
table(quine$c.Days)


# split quine in two subsets
suppressWarnings(RNGversion("3.5.0"))
set.seed(1111)
lab.A <- sample(nrow(quine), 70, replace=TRUE)
quine.A <- quine[lab.A, 1:4]
quine.B <- quine[-lab.A, c(1:3,6)]

# compute the tables required by Fbwidths.by.x()
freq.xA <- xtabs(~Eth+Sex+Age, data=quine.A)
freq.xB <- xtabs(~Eth+Sex+Age, data=quine.B)

freq.xy <- xtabs(~Eth+Sex+Age+Lrn, data=quine.A)
freq.xz <- xtabs(~Eth+Sex+Age+c.Days, data=quine.B)

# apply Fbwidths.by.x()
bb <- Fbwidths.by.x(tab.x=freq.xA+freq.xB, 
                           tab.xy=freq.xy,  tab.xz=freq.xz,
                           warn=FALSE)
bb$sum.unc
cc <- selMtc.by.unc(tab.x=freq.xA+freq.xB, 
                           tab.xy=freq.xy,  tab.xz=freq.xz, corr.d=0)
cc$ini.ord
cc$av.df


}

\keyword{multivariate}