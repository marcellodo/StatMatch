% --- Source file:  ---
\name{Frechet.bounds.cat}
\alias{Frechet.bounds.cat}
\title{Frechet bounds of cells in a contingency table}

\description{
This function permits to derive the bounds for cell probabilities of the table Y vs. Z starting from the marginal tables (\bold{X} vs. Y), (\bold{X} vs. Z) and the joint distribution of the \bold{X} variables.
}

\usage{
Frechet.bounds.cat(tab.x, tab.xy, tab.xz, print.f="tables", tol= 0.0001) 
}

\arguments{

\item{tab.x}{
A \R table crossing the \bold{X} variables.  This table must be obtained by using the function \code{\link[stats]{xtabs}} or \code{\link[base]{table}},  e.g. \cr
\code{tab.x <- xtabs(~x1+x2+x3, data=data.all)}. \cr
When \code{tab.x = NULL} then \code{tab.xy} and \code{tab.xz} must be supplied.
}

\item{tab.xy}{
A \R table of \bold{X} vs. Y variable.  This table must be obtained by using the function \code{\link[stats]{xtabs}} or \code{\link[base]{table}}, e.g. \cr
\code{table.xy <- xtabs(~x1+x2+x3+y, data=data.A)}.

A single categorical Y variable is allowed.  One or more categorical variables can be considered as \bold{X} variables (common variables).  Obviously, the same \bold{X} variables in \code{tab.x} must be available in \code{tab.xy}.  Moreover, it is assumed that the joint distribution of the \bold{X} variables computed from \code{tab.xy} is equal to \code{tab.x}; a warning appears if this is not true (see argument \code{tol}).

When \code{tab.x = NULL} then \code{tab.xy} should be a one--dimensional table providing the marginal distribution of the Y variable.
}

\item{tab.xz}{
A \R table of \bold{X} vs. Z variable.  This table must be obtained by using the function \code{\link[stats]{xtabs}} or \code{\link[base]{table}}, e.g. \cr
\code{tab.xz <- xtabs(~x1+x2+x3+z, data=data.B)}.

A single categorical Z variable is allowed.  One or more categorical variables can be considered as \bold{X} variables (common variables).  The same \bold{X} variables in \code{tab.x} must be available in \code{tab.xz}.  Moreover, it is assumed that the joint distribution of the \bold{X} variables computed from \code{tab.xz} is equal to \code{tab.x}; a warning appears if this is not true (see argument \code{tol}).

When \code{tab.x = NULL} then \code{tab.xz} should be a one--dimensional table providing the marginal distribution of the Z variable.

}

\item{print.f}{
A string: when \code{print.f="tables"} (default) all the cells' estimates will be saved as tables in a list.  On the contrary, if \code{print.f="data.frame"}, they will be saved as columns of a data.frame.
}

\item{tol}{
Tolerance used in comparing joint distributions as far as \bold{X} variables are considered (default \code{tol= 0.0001}); the joint distribution of the \bold{X} variables computed from \code{tab.xy} and \code{tab.xz} should be equal to that in \code{tab.x}.
}

}

\details{
This function permits to compute the Frechet bounds for the relative frequencies in the contingency table of Y vs. Z, starting from the distributions P(Y|X), P(Z|X) and P(X).  The bounds for the relative frequencies \eqn{p_{j,k}}{p(y=j,z=k)} in the table Y vs. Z are:

\deqn{  p^{(low)}_{Y=j,Z=k} = \sum_{i}  p_{X=i}\max (0; p_{Y=j|X=i} + p_{Z=k|X=i} - 1 ) }{p(Y=j,Z=k) >= sum_i(p(X=i) * max(0; p(Y=j|X=i) + p(Z=k|X=i) - 1) )}

\deqn{ p^{(up)}_{Y=j,Z=k} = \sum_{i} p_{X=i} \min ( p_{Y=j|X=i};  p_{Z=k|X=i})}{p(Y=j,Z=k) <= sum_i(p(X=i) * min(p(Y=j|X=i),p(Z=k|X=i))) }


The relative frequencies \eqn{p_{X=i}=n_i/n}{p(X=i)=n_i/n} are computed from the frequencies in \code{tab.x}; \cr
the relative frequencies \eqn{p_{Y=j|X=i}=n_{ij}/n_{i+}}{p(Y=j|X=i)=n_ij/n_i.} are computed from the \code{tab.xy}, \cr
finally, \eqn{p_{Z=k|X=i}=n_{ik}/n_{k+}}{p(Z=k|X=i)=n_ik/n_i.} are derived from \code{tab.xy}.

It is assumed that the marginal distribution of the \bold{X} variables is the same in all the input tables: \code{tab.x}, \code{tab.xy} and \code{tab.xz}.  If this is not true a warning message will appear.

Note that the cells bounds for the relative frequencies in the contingency table of Y vs. Z are computed also without considering the \bold{X} variables:

\deqn{ \max\{0; p_{Y=j} + p_{Z=k} - 1\} \leq p_{Y=j,Z=k} \leq \min \{ p_{Y=j}; p_{Z=k}\}}{max(0;p(Y=j)+p(Z=k)-1) <= p(Y=j,Z=k) <= min(p(Y=j);p(Z=k))}

These bounds represent the unique output when \code{tab.x = NULL}.

Finally, the contingency table of Y vs. Z estimated under the Conditional Independence Assumption (CIA) is obtained by considering:

\deqn{ p_{Y=j,Z=k} = p_{Y=j|X=i} \times p_{Z=k|X=i} \times p_{X=i}.}{p(Y=i,Z=k) = p(Y=j|X=i)*p(Z=k|X=i)*p(X=i)}

When \code{tab.x = NULL} then it is also provided the expected table under the assumption of independence between Y and Z:

\deqn{ p_{Y=j,Z=k} = p_{Y=j} \times p_{Z=k}.}{p(Y=i,Z=k) = p(Y=j)*p(Z=k)*}

}

\value{
When \code{print.f="tables"} (default) a list with the following components:

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

\item{uncertainty}{
The uncertainty associated to input data, summarized in terms average width of uncertainty bounds with and without conditioning on the \bold{X} variables estimated, overall uncertainty estimated according to the suggestion in Conti \emph{et al.} (2012) (see \code{\link[StatMatch]{Fbwidths.by.x}} for major details).
}


When \code{print.f="data.frame"} the output list contains just two components: 

\item{bounds}{
A data.frame whose columns reports the estimated uncertainty bounds.
}

\item{uncertainty}{
The uncertainty associated to input data, summarized in terms average width of uncertainty bounds with and without conditioning on the \bold{X} variables estimated, overall uncertainty estimated according to the suggestion in Conti \emph{et al.} (2012) (see \code{\link[StatMatch]{Fbwidths.by.x}} for major details).
}

}

\references{

Ballin, M., D'Orazio, M., Di Zio, M., Scanu, M. and Torelli, N. (2009) \dQuote{Statistical Matching of Two Surveys with a Common Subset}. \emph{Working Paper}, \bold{124}. Dip. Scienze Economiche e Statistiche, Univ. di Trieste, Trieste. 

Conti P.L, Marella, D., Scanu, M. (2012) \dQuote{Uncertainty Analysis in Statistical Matching}. \emph{Journal of Official Statistics}, \bold{28}, pp. 69--88.

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
set.seed(7654)
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

#compare marg. distribution of Xs in A and B
comp.prop(p1=margin.table(freq.xy,c(1,2)), p2=margin.table(freq.xz,c(1,2)), 
          n1=nrow(quine.A), n2=nrow(quine.B))

# harmonize distr. of Sex vs. Age before applying
# Frechet.bounds.cat()

N <- nrow(quine)
quine.A$pop <- N
quine.A$f <- N/70 # reciprocal sampling fraction
quine.B$pop <- N
quine.B$f <- N/(N-70)

# derive the table of Sex vs. Age related to the whole data set
tot.sex.age <- colSums(model.matrix(~Sex*Age-1, data=quine))
tot.sex.age

# use hamonize.x() to harmonize the Sex vs. Age between
# quine.A and quine.B

# create svydesign objects
require(survey)
svy.qA <- svydesign(~1, weights=~f, fpc=~pop, data=quine.A)
svy.qB <- svydesign(~1, weights=~f, fpc=~pop, data=quine.B)

# apply harmonize.x 
out.hz <- harmonize.x(svy.A=svy.qA, svy.B=svy.qB, form.x=~Sex*Age-1, x.tot=tot.sex.age)

# compute the new tables required by Frechet.bounds.cat()
freq.x <- xtabs(out.hz$weights.A~Sex+Age, data=quine.A)
freq.xy <- xtabs(out.hz$weights.A~Sex+Age+Eth, data=quine.A)
freq.xz <- xtabs(out.hz$weights.B~Sex+Age+Lrn, data=quine.B)

#compare marg. distribution of Xs in A and B
comp.prop(p1=margin.table(freq.xy,c(1,2)), p2=margin.table(freq.xz,c(1,2)), 
          n1=nrow(quine.A), n2=nrow(quine.B))

# apply Frechet.bounds.cat()
bounds.yz <- Frechet.bounds.cat(tab.x=freq.x, tab.xy=freq.xy,
        tab.xz=freq.xz, print.f="data.frame")
bounds.yz

}


\keyword{multivariate}