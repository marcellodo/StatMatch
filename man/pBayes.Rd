\name{pBayes}
\alias{pBayes}
\title{Pseudo-Bayes estimates of cell probabilities}

\description{
  Estimation of cells counts using the pseudo-Bayes estimator.
}

\usage{
pBayes(x, method="m.ind", const=NULL)
}

\arguments{

\item{x}{
 A contingency table with observed cell counts. Typically the output of \code{\link[base]{table}} or \code{\link[stats]{xtabs}}. More in general an R \code{array} with the counts.
}

\item{method}{ 
 The method for estimating the final cell frequencies. The following options are available:

\code{method = "Jeffreys"}, consists in adding 0.5 to each cell before estimation of the relative frequencies. 

\code{method = "minimax"}, consists in adding \eqn{\sqrt(n)/c}{sqrt(n)/c} to each cell before estimation of the relative frequencies, being \eqn{n}{n} the sum of all the counts and \eqn{c}{c} the number of cells in the table.

\code{method = "invcat"}, consists in adding \eqn{1/c}{1/c} to each cell before estimation of the relative frequencies.

\code{method = "user"}, consists in adding a used defined constant \eqn{a}{a} (\eqn{a>0}{a>0}) to each cell before estimation of the relative frequencies. The constant \eqn{a}{a} should be passed via the argument \code{const}.

\code{method = "m.ind"}, the prior guess for the unknown cell probabilities is obtained by considering estimated probabilities under the mutual independence hypothesis. This option is available when dealing with at least two-way contingency tables \cr (\code{length(dim(x))>=2}).

\code{method = "h.assoc"}, the prior guess for the unknown cell probabilities is obtained by considering estimated probabilities under the homogeneous association hypothesis. This option is available when dealing with at least two-way contingency tables (\code{length(dim(x))>=2}).

}

\item{const}{
 Numeric value, a user defined constant \eqn{a}{a} (\eqn{a>0}{a>0}) added to each cell before estimation of the relative frequencies when \code{method = "user"}. As a general rule of thumb, it is preferable to avoid that the sum of constant over all the cells in table is greater than \eqn{0.20 \times n}{0.20*n}.
}

}

\details{
This function estimates the frequencies in a contingency table by using the pseudo-Bayes approach. In practice the estimator being considered is a weighted average of the observed cells counts \eqn{n_h}{n_h} and a suitable prior guess for cells probabilities \eqn{\gamma_h}{gamma_h}:

\deqn{\tilde{p}_h = \frac{n}{n+K} \hat{p}_h + \frac{K}{n+K} \gamma_h }{ep_h = n/(n+K)*p_h + K/(n+K)*gamma_h }

\eqn{K}{K} depends on the parameters of Dirichlet prior distribution being considered (for major details see Chapter 12 in Bishop et al., 1974).
It is worth noting that with a constant prior guess \eqn{\gamma_h=1/c}{gamma_h=1/c} (\eqn{h=1,2,\cdots, c}{h=1,2,...,c}), then \eqn{K=1}{K=1} and in practice corresponds to adding \eqn{1/c}{1/c} to each cell before estimation of the relative frequencies (\code{method = "invcat"}); \eqn{K=c/2}{K=c/2} when the constant 0.5 is added to each cell (\code{method = "Jeffreys"}); finally \eqn{K=\sqrt{n}}{sqrt(n)} when the quantity \eqn{\sqrt{n}/c}{sqrt(n)/c} is added to each cell (\code{method = "minimax"}). All these cases corresponds to adding a flattening constant; the higher is the value of \eqn{K}{K} the more the estimates will be shrinked towards \eqn{\gamma_h=1/c}{gamma_h=1/c} (flattening).

When \code{method = "m.ind"} the prior guess \eqn{\gamma_h}{gamma_h} is estimated under the hypothesis of mutual independence between the variables crossed in the initial contingency table \code{x}, supposed to be at least a two-way table.  In this case the value of \eqn{K}{K} is estimated via a data driven approach by considering

\deqn{ \hat{K} = \frac{1 - \sum_{h} \hat{p}_h^2}{\sum_{h} \left( \hat{\gamma}_h - \hat{p}_h \right)^2 } }{ eK = (1 - sum(p_h^2))/(sum(gamma_h - p_h)^2) }

On the contrary, when \code{method = "h.assoc"} the prior guess \eqn{\gamma_h}{gamma_h} is estimated under the hypothesis of homogeneous association between the variables crossed in the initial contingency table \code{x}. 

Please note that when the input table is estimated from sample data where a weight is assigned to each unit, the weights should be used in estimating the input table, but it is suggested to rescale them so that they sum up to \emph{n}, the sample size.
  
}

\value{
A \code{list} object with three components. 

\item{info}{A vector with the sample size \code{"n"}, the number of cells (\code{"no.cells"}) in \code{x}, the average cell frequency (\code{"av.cfr"}), the number of cells showing frequencies equal to zero (\code{"no.0s"}), the \code{const} input argument, the chosen/estimated \eqn{K}{K} (\code{"K"}) and the relative size of \eqn{K}{K}, i.e. \eqn{K/(n+K)}{K/(n+K)} (\code{"rel.K"}).}

\item{prior}{ A table having the same dimension as \code{x} with the considered prior values for the cell frequencies.}

\item{pseudoB}{ A table with having the same dimension as \code{x} providing the pseudo-Bayes estimates for the cell frequencies in \code{x}.}


}

\references{

Bishop Y.M.M., Fienberg, S.E., Holland, P.W. (1974) \emph{Discrete Multivariate Analysis: Theory and Practice.} The Massachusetts Institute of Technology

}

\author{
 Marcello D'Orazio \email{mdo.statmatch@gmail.com} 
}

% \seealso{ 
%
% }

\examples{

data(samp.A, package="StatMatch")
tab <- xtabs(~ area5 + urb + c.age + sex + edu7, data = samp.A)
out.pb <- pBayes(x=tab, method="m.ind")
out.pb$info

out.pb <- pBayes(x=tab, method="h.assoc")
out.pb$info

out.pb <- pBayes(x=tab, method="Jeffreys")
out.pb$info

# usage of weights in estimating the input table
n <- nrow(samp.A)
r.w <- samp.A$ww / sum(samp.A$ww) * n   # rescale weights to sum up to n 
tab.w <- xtabs(r.w ~ area5 + urb + c.age + sex + edu7, data = samp.A)
out.pbw <- pBayes(x=tab.w, method="m.ind")
out.pbw$info

}

\keyword{multivariate}
