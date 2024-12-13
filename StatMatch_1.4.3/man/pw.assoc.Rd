\name{pw.assoc}
\alias{pw.assoc}
\title{Pairwise measures between categorical variables}

\description{
  This function computes some association and \emph{Proportional Reduction in Error} (PRE) measures between a categorical nominal variable and each of the other available predictors (being also categorical variables).
}

\usage{
pw.assoc(formula, data, weights=NULL, out.df=FALSE)
}

\arguments{

\item{formula}{
 A formula of the type \code{y~x1+x2} where \code{y} denotes the name of the categorical variable (a \code{factor} in \pkg{R}) which plays the role of the dependent variable, while \code{x1} and \code{x2} are the name of the predictors (both categorical variables).  Numeric variables are not allowed; eventual numerical variables should be categorized (see function \code{\link[base]{cut}}) before being passed to \code{pw.assoc}.
}

\item{data}{ 
The data frame which contains the variables called by \code{formula}.
}

\item{weights}{
 The name of the variable in \code{data} which provides the units' weights.  Weights are used to estimate frequencies (a cell frequency is estimated by summing the weights of the units which present the given characteristic).  Default is \code{NULL} (no weights available) and each unit counts 1. When case weight are provided, then they are scales so that their sum equals \emph{n}, the sample size (assumed to be \code{nrow(data)}).
}

\item{out.df}{
 Logical. If \code{NULL} measures will be organized in a data frame (a column for each measure).
}

}

\details{
This function computes some association, PRE measures, AIC and BIC for each couple response-predictor that can be created  starting from argument \code{formula}. In particular, a two-way contingency table \eqn{X \times Y}{X x Y} is built for each available X variable (X in rows and Y in columns); then the following measures are considered.

Cramer's \emph{V}: 

\deqn{ V=\sqrt{\frac{\chi^2}{n \times min\left[I-1,J-1\right]} } }{ (Chi^2/(N*min(I-1,J-1)))^0.5}

\emph{n} is the sample size, \emph{I} is the number of rows (categories of X) and \emph{J} is the number of columns (categories of Y).  Cramer's \emph{V} ranges from 0 to 1. 

Bias-corrected Cramer's \emph{V} (\eqn{V_c}{V_c}) proposed by Bergsma (2013).

Mutual information:

\deqn{ I(X;Y) = \sum_{i,j} p_{ij} \, log \left( \frac{p_{ij}}{p_{i+} p_{+j}} \right) }{ sum_ij (p_ij*log p_ij/(p_i+*p_+j) ) }  

equal to 0 in case of independence but with infinite upper bound, i.e.  \eqn{0 \leq I(X;Y) < \infty}{0 <= I(X;Y) < Infinite}. In it \eqn{p_{ij}=n_{ij}/n }{p_ij=n_ij/n}.

% A normalized version of \emph{I} also known as (symmetric) Theil's \emph{U} is:
%
% \deqn{ U = 2\frac{I}{- \sum_i p_{i+} log p_{i+} - \sum_j p_{+j} log p_{+j}  }}{ U = 2*I / ( H_Y + H_X) }
%
% and ranges from 0 (independence) to 1. 

A normalized version of \eqn{I(X;Y)}{I(X;Y)}, ranging from 0 (independence) to 1 and not affected by number of categories (\emph{I} and \emph{J}):

\deqn{I(X;Y)^* = \frac{I(X;Y)}{min(H_X, H_Y) } }{I* = I/(min(H_X, H_Y) )}

being \eqn{H_X}{H_X} and \eqn{H_Y}{H_Y} the entropy of the variable X and Y, respectively.


Goodman-Kruskal \eqn{\lambda(Y|X)}{lambda(Y|X)} (i.e. response conditional on the given predictor):

\deqn{ \lambda(Y|X) = \frac{\sum_{i=1}^I max_{j}(p_{ij}) - max_{j}(p_{+j})}{1-max_{j}(p_{+j})} }{ lambda(Y|X) = (sum_i max_j(p_ij) - max_j(p_+j))/(1 - max_j(p_+j)) }

It ranges from 0 to 1, and denotes how much the knowledge of the row variable X (predictor) helps in reducing the prediction error of the values of the column variable Y (response).


Goodman-Kruskal \eqn{\tau(Y|X)}{tau(Y|X)}:

\deqn{ \tau(Y|X) = \frac{ \sum_{i=1}^I \sum_{j=1}^J p^2_{ij}/p_{i+} - \sum_{j=1}^J p_{+j}^2}{1 -  \sum_{j=1}^J p_{+j}^2} }{ tau(Y|X) = (sum_ij  p^2_ij / p_i+ - sum_j p^2_+j)/(1 - sum_j p^2_+j)  }

It takes values in the interval [0,1] and has the same PRE meaning of the lambda.

Theil's uncertainty coefficient:

\deqn{ U(Y|X) = \frac{\sum_{i=1}^I \sum_{j=1}^J p_{ij} log(p_{ij}/p_{i+}) - \sum_{j=1}^J p_{+j} log p_{+j}}{- \sum_{j=1}^J p_{+j} log p_{+j}} }{ U(Y|X) = (sum_ij p_ij log (p_ij/pi+) - sum_j p_+j log p_+j) / (- sum_j p_+j log p_+j) }

It takes values in the interval [0,1] and measures the reduction of uncertainty in the column variable Y due to knowing the row variable X. Note that the numerator of U(Y|X) is the mutual information I(X;Y)

It is worth noting that \eqn{\lambda}{lambda}, \eqn{\tau}{tau} and \emph{U} can be viewed as measures of the proportional reduction of the variance of the Y variable when passing from its marginal distribution to its conditional distribution given the predictor X, derived from the general expression (cf. Agresti, 2002, p. 56):

\deqn{ \frac{V(Y) - E[V(Y|X)]}{V(Y)}}{ (V(Y) - E[V(Y|X)])/V(Y) }

They differ in the way of measuring variance, in fact it does not exist a general accepted definition of the variance for a categorical variable.

Finally, AIC (and BIC) is calculated, as suggested in Sakamoto and Akaike (1977). In particular:

\deqn{ AIC(Y|X) = -2 \sum_{i,j} n_{ij} \, log \left( \frac{n_{ij}}{n_{i+}} \right) + 2I(J - 1) }{-2*sum_ij (n_ij * log(n_ij/n_i+) +2*I*(J-1) )  }

\deqn{ BIC(Y|X) = -2 \sum_{i,j} n_{ij} \, log \left( \frac{n_{ij}}{n_{i+}} \right) +I(J-1) log(n)    }{-2*sum_ij (n_ij * log(n_ij/n_i+) +I*(J-1)*log(n) )  }

being \eqn{I(J-1)}{I*(J-1)} the parameters (conditional probabilities) to estimate. Note that the \R package \pkg{catdap} provides functions to identify the best subset of predictors based on AIC.

Please note that the missing values are excluded from the tables and therefore excluded from the estimation of the various measures.
}

\value{
When \code{out.df=FALSE} (default) a \code{list} object with four components: 

\item{V}{A vector with the estimated Cramer's V for each couple response-predictor.}

\item{bcV}{A vector with the estimated bias-corrected Cramer's V for each couple response-predictor.}

\item{mi}{A vector with the estimated mutual information I(X;Y) for each couple response-predictor.}

\item{norm.mi}{A vector with the normalized mutual information I(X;Y)* for each couple response-predictor.}

\item{lambda}{ A vector with the values of Goodman-Kruscal \eqn{\lambda(Y|X)}{lambda(Y|X)} for each couple response-predictor.}

\item{tau}{ A vector with the values of Goodman-Kruscal \eqn{\tau(Y|X)}{tau(Y|X)} for each couple response-predictor.}

\item{U}{ A vector with the values of Theil's uncertainty coefficient U(Y|X) for each couple response-predictor.}

\item{AIC}{ A vector with the values of AIC(Y|X) for each couple response-predictor.}

\item{BIC}{ A vector with the values of BIC(Y|X) for each couple response-predictor.}

\item{npar}{ A vector with the number of parameters (conditional probabilities) estimated to calculate AIC and BIC for each couple response-predictor.}

When \code{out.df=TRUE} the output will be a data.frame with a column for each measure. 

}

\references{

Agresti A (2002) \emph{Categorical Data Analysis. Second Edition}. Wiley, new York.

Bergsma W (2013) A bias-correction for Cramer's V and Tschuprow's T. \emph{Journal of the Korean Statistical Society}, 42, 323--328.

The Institute of Statistical Mathematics (2018). catdap: Categorical Data  Analysis Program Package. R package version 1.3.4. \url{https://CRAN.R-project.org/package=catdap}

Sakamoto Y and Akaike, H (1977) Analysis of Cross-Classified Data by AIC. \emph{Ann. Inst. Statist. Math.}, 30, 185-197.
}

\author{
 Marcello D'Orazio \email{mdo.statmatch@gmail.com} 
}

% \seealso{ 
%
% }

\examples{
data(quine, package="MASS") #loads quine from MASS
str(quine)

# how Lrn is response variable
pw.assoc(Lrn~Age+Sex+Eth, data=quine)

# usage of units' weights
quine$ww <- runif(nrow(quine), 1,4) #random gen  1<=weights<=4
pw.assoc(Lrn~Age+Sex+Eth, data=quine, weights="ww")

}

\keyword{multivariate}
