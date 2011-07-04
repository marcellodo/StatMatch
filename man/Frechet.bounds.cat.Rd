% --- Source file:  ---
\name{Frechet.bounds.cat}
\alias{Frechet.bounds.cat}
\title{Frechet bounds of cells in a contingency table}

\description{
This function permits to derive the bounds for cell probabilities of the table Y vs. Z starting from the marginal tables (\bold{X} vs. Y), (\bold{X} vs. Z) and the joint distribution of the \bold{X} variables.
}

\usage{
Frechet.bounds.cat(tab.x, tab.xy, tab.xz, print.f="tables") 
}

\arguments{

\item{tab.x}{
A \R table crossing the \bold{X} variables.  This table must be obtained by using the function \code{\link[stats]{xtabs}} or \code{\link[base]{table}},  e.g. \cr
\code{tab.x <- xtabs(~x1+x2+x3, data=data.all)}.
}

\item{tab.xy}{
A \R table of \bold{X} vs. Y variable.  This table must be obtained by using the function \code{\link[stats]{xtabs}} or \code{\link[base]{table}}, e.g. \cr
\code{table.xy <- xtabs(~x1+x2+x3+y, data=data.A)}.

A single categorical Y variable is allowed.  One or more categorical variables can be considered as \bold{X} variables (common variables).  Obviously, the same \bold{X} variables in \code{tab.x} must be available in \code{tab.xy}.  Moreover, it is assumed that the joint distribution of the \bold{X} variables computed from \code{tab.xy} is equal to \code{tab.x}; a warning appears if this is not true.
}

\item{tab.xz}{
A \R table of \bold{X} vs. Z variable.  This table must be obtained by using the function \code{\link[stats]{xtabs}} or \code{\link[base]{table}}, e.g. \cr
\code{tab.xz <- xtabs(~x1+x2+x3+z, data=data.B)}.

A single categorical Z variable is allowed.  One or more categorical variables can be considered as \bold{X} variables (common variables).  The same \bold{X} variables in \code{tab.x} must be available in \code{tab.xz}.  Moreover, it is assumed that the joint distribution of the \bold{X} variables computed from \code{tab.xz} is equal to \code{tab.x}; a warning appears if this is not true.
}

\item{print.f}{
A string specifying the data structure of the output. When \cr \code{print.f="tables"} (default) all the results will be saved as tables in a list.  On the contrary, if \code{print.f="data.frame"}, all results will be saved as columns of a data.frame.
}

}

\details{
This function permits to compute the Frechet bounds for the relative frequencies in the contingency table of Y vs.Z, starting from the distributions P(Y|X), P(Z|X) and P(X).  The bounds for the relative frequencies \eqn{p_{j,k}}{p(y=j,z=k)} in the table Y vs. Z are:

\deqn{  p^{low}_{YZ}(j,k) = \sum_{i}  p_X(i)\max (0; p_{Y|X}(j|i) + p_{Z|X}(k|i)-1 ) }{p(y=j,z=k) >= sum_i(p(x=i) * max(0; p(y=j|x=i) + p(z=k|x=i) - 1) )}

\deqn{ p^{up}_{YZ}(j,k) = \sum_{i} p_X(i) \min ( p_{Y|X}(j|i);  p_{Z|X}(k|i))}{p(y=j,z=k) <= sum_i(p(x=i) * min(p(y=j|x=i),p(z=k|x=i))) }


The relative frequencies \eqn{p_X(i)=n_i/n}{p(x=i)=n_i/n \bullet} are computed from the frequencies in \code{tab.x}; \cr
the relative frequencies \eqn{p_{Y|X}(j|i)=n_{ij}/n_{i \bullet}}{p(y=j|x=i)=n_ij/n_i.} are computed from the \code{tab.xy}, \cr
finally, \eqn{p_{Z|X}(k|i)=n_{ik}/n_{k \bullet}}{p(z=k|x=i)=n_ik/n_i.} are derived from \code{tab.xy}.

It is assumed that the marginal distribution of the \bold{X} variables is the same in all the input tables: \code{tab.x}, \code{tab.xy} and \code{tab.xz}.  If this is not true a warning message will appear.

Note that the cells bounds for the relative frequencies in the contingency table of Y vs. Z are computed also without considering the \bold{X} variables:

\deqn{ \max\{0; p_{Y}(j) + p_{Z}(k)-1\} \leq p_{YZ}(j,k) \leq \min \{ p_{Y}(j); p_{Z}(k)\}}{max(0;p(y=j)+p(z=k)-1) <= p(y=j,z=k) <= min(p(y=j);p(z=k))}

Finally, the contingency table of Y vs. Z estimated under the Conditional Independence Assumption (CIA) is obtained by considering:

\deqn{ p_{YZ}(j,k) = p_{Y|X}(j|i) \times p_{Z|X}(k|i) \times p_{X}(i).}{p(y=i,z=k) = p(y=j|x=i)*p(z=k|x=i)*p(x=i)}

}

\value{
When \code{print.f="tables"} (default) a list with the following tables:

\item{low.u}{
The estimated lower bounds for the relative frequencies in the table Y vs. Z without conditioning on the \bold{X} variables.
}

\item{up.u}{
The estimated upper bounds for the relative frequencies in the table Y vs. Z without conditioning on the \bold{X} variables.
}

\item{CIA}{
The estimated relative frequencies in the table Y vs. Z under the Conditional Independence Assumption (CIA).
}

\item{low.cx}{
The estimated lower bounds for the relative frequencies in the table Y vs. Z when conditioning on the \bold{X} variables.
}

\item{up.cx}{
The estimated upper bounds for the relative frequencies in the table Y vs. Z when conditioning on the \bold{X} variables.
}

When \code{print.f="data.frame"} the estimated tables are saved as columns of a data.frame.

}

\references{

Ballin, M., D'Orazio, M., Di Zio, M., Scanu, M. and Torelli, N. (2009) \dQuote{Statistical Matching of Two Surveys with a Common Subset}. \emph{Working Paper}, \bold{124}. Dip. Scienze Economiche e Statistiche, Univ. di Trieste, Trieste. 

D'Orazio, M., Di Zio, M. and Scanu, M. (2006). \emph{Statistical Matching: Theory and Practice.} Wiley, Chichester.

}


\author{
 Marcello D'Orazio \email{madorazi@istat.it} 
}

\seealso{ 
\code{\link[StatMatch]{Fbwidths.by.x}}, \code{\link[StatMatch]{harmonize.x}}
}

\examples{
data(quine, package="MASS") #loads quine from MASS
str(quine)

# split quine in two subsets
set.seed(765)
lab.A <- sample(nrow(quine), 70, replace=TRUE)
quine.A <- quine[lab.A, 1:3]
quine.B <- quine[-lab.A, 2:4]

# compute the tables required by Frechet.bounds.cat()
freq.x <- xtabs(~Sex+Age, data=quine.A)
freq.xy <- xtabs(~Sex+Age+Eth, data=quine.A)
freq.xz <- xtabs(~Sex+Age+Lrn, data=quine.B)

# apply Frechet.bounds.cat()
bounds.yz <- Frechet.bounds.cat(tab.x=freq.x, tab.xy=freq.xy,
        tab.xz=freq.xz, print.f="data.frame")
bounds.yz

# harmonize distr. of Sex vs. Age before applying
# Frechet.bounds.cat()

quine.A$f <- 70/nrow(quine) # sampling fraction
quine.B$f <- (nrow(quine)-70)/nrow(quine)

# derive the table of Sex vs. Age related to the whole data set
tot.sex.age <- xtabs(~Sex+Age, data=quine)
tot.sex.age

# use hamonize.x() to harmonize the Sex vs. Age between
# quine.A and quine.B

# create svydesign objects
svy.qA <- svydesign(~1, fpc=~f, data=quine.A)
svy.qB <- svydesign(~1, fpc=~f, data=quine.B)

# apply harmonize.x using poststratification
out.hz <- harmonize.x(svy.A=svy.qA, svy.B=svy.qB, form.x=~Sex+Age,
          cal.method="poststratify")

# compute the new tables required by Frechet.bounds.cat()
freq.x <- xtabs(out.hz$weights.A~Sex+Age, data=quine.A)
freq.xy <- xtabs(out.hz$weights.A~Sex+Age+Eth, data=quine.A)
freq.xz <- xtabs(out.hz$weights.B~Sex+Age+Lrn, data=quine.B)

# apply Frechet.bounds.cat()
bounds.yz <- Frechet.bounds.cat(tab.x=freq.x, tab.xy=freq.xy,
        tab.xz=freq.xz, print.f="data.frame")
bounds.yz

}


\keyword{multivariate}