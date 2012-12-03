\name{pw.assoc}
\alias{pw.assoc}
\title{Pairwise association measure between categorical variables}

\description{
  This function computes some association measures between a categorical nominal variable and each of the other  available predictors (also categorical variables).
}

\usage{
pw.assoc(formula, data, weights=NULL, freq0c=NULL)
}

\arguments{

\item{formula}{
 A formula of the type \code{y~x1+x2} where \code{y} denotes the name of the categorical variable (a \code{factor} in \pkg{R}) which plays the role of the dependent variable while \code{x1} and \code{x2} are the name of the predictors (both categorical variables).  Numeric variables are not allowed; eventual numerical variables should be categorized (see function \code{\link[base]{cut}}) before being passed to \code{pw.assoc}.
}

\item{data}{ 
The data frame which contains the variables called by \code{formula}.
}

\item{weights}{
 The name of the eventual variable in \code{data} which provides the units' weights.  Weights are used to estimate frequencies (a cell frequency is estimated by summing the weights of the units which present the given characteristics).  Default is \code{NULL} (no weights available, each unit counts 1).
}

\item{freq0c}{
 A small number which is substituted to eventual cells with zero frequencies in order to avoid computation failures.  When \code{NULL} (default) a cell with zero frequency is substitutes with 1/N^2, being N the sample size.
}

}

\details{
This function computes some association measures among the response variable and each of the predictors specified in the \code{formula}. The following association measure are considered:

Cramer's \emph{V}: 

\deqn{ V=\sqrt{\frac{\chi^2}{N \times min\left[I-1,J-1\right]} } }{ (Chi^2/(N*min(I-1,J-1)))^0.5}

N is the sample size, I is the number of rows and J is the number of columns.  Cramer's V ranges from 0 to 1. 

Goodman--Kruskal \eqn{\lambda(R|C)}{lambda(R|C)}:

\deqn{ \lambda(R|C) = \frac{\sum_{j=1}^J max_{i}(p_{ij}) - max_{i}(p_{i+})}{1-max_{i}(p_{i+})} }{ lambda(R|C) = (sum_j max_i(p_ij) - max_i(p_i+))/(1 - max_i(p_i+)) }

It ranges from 0 to 1, and denotes how much the knowledge of the column variable (predictor) helps in reducing the prediction error of the values of the row variable.


Goodman--Kruskal \eqn{\tau(R|C)}{tau(R|C)}:

\deqn{ \tau(R|C) = \frac{ \sum_{i=1}^I \sum_{j=1}^J p^2_{ij}/p_{+j} - \sum_{i=1}^I p_{i+}^2}{1 -  \sum_{i=1}^I p_{i+}^2} }{ tau(R|C) = (sum_ij  p^2_ij / p_+j - sum_i p_i+)/(1 - sum_i p_i+)  }

It takes values in the interval [0,1] and has the same PRE meaning of the lambda.

Theil's Uncertainty coefficient:

\deqn{ U(R|C) = \frac{\sum_{i=1}^I \sum_{j=1}^J p_{ij} log(p_{ij}/p_{+j}) - \sum_{i=1}^I p_{i+} log p_{i+}}{- \sum_{i=1}^I p_{i+} log p_{i+}} }{ U(R|C) = (sum_ij p_ij log (p_ij/p+j) - sum_i p_i+ log p_i+) / (-sum_i p_i+ log p_i+) }

It takes values in the interval [0,1] and measure the reduction of uncertainty in the row variable due to knowing the column variable.


It is worth noting that \eqn{\lambda}{lambda}, \eqn{\tau}{lambda} and U are asymmetric measures of the proportional reduction of the variance of the row column when passing from its marginal distribution to its conditional distribution given the column variable obtained  starting from the general expression (cf. Agresti, 2002, p. 56):

\deqn{ \frac{V(R) - E[V(R|C)]}{V(R)}}{ (V(R) - E[V(R|C)])/V(R) }

They differ in the way of measuring variance, in fact it does not exist a general accepted definition of the variance of a categorical variable.

}

\value{
A \code{list} object with for components. 

\item{V}{A vector with the estimated Cramer's V for each couple response-predictor.}

\item{labda}{ A vector with the values of Goodman-Kruscal \eqn{\lambda(R|C)}{lambda(R|C)} for each couple response-predictor.}

\item{tau}{ A vector with the values of Goodman-Kruscal \eqn{\tau(R|C)}{tau(R|C)} for each couple response-predictor.}

\item{U}{ A vector whit the values of Theil's uncertainty coefficient U(R|C) for each couple response-predictor.}

}

\references{

Agresti A (2002) \emph{Categorical Data Analysis. Second Edition}. Wiley, new York.

}

\author{
 Marcello D'Orazio \email{madorazi@istat.it} 
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
