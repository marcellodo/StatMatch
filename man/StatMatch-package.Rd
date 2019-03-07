\name{StatMatch-package}
\alias{StatMatch-package}
\alias{StatMatch}
\docType{package}
\title{
Statistical Matching or Data Fusion
}
\description{
Functions to perform statistical matching (aka data fusion), i.e. the integration of two data sources.  Some functions can also be used to impute missing values in data sets through hot deck imputation methods.

}
\details{

\emph{Statistical matching} (aka \emph{data fusion}) aims at integrating two data sources (typically samples)  referred to same target population and sharing a number of variables.  The final objective is that of studying the relationship of variables not jointly observed in a single data sources.  The integration can be performed at \emph{micro} (a \emph{synthetic} or \emph{fused} file is the output) or \emph{macro} level (goal is estimation of correlation coefficients, regression coefficients, contingency tables, etc.).  
Nonparametric hot deck imputation methods (random, rank and nearest neighbour donor) can be used to derive the synthetic data set.  In alternative it is possible to use a mixed (parametric-nonparametric) procedure based on predictive mean matching.  Methods to perform statistical matching when dealing with data from complex sample surveys are available too.  Finally some functions can be used to explore the uncertainty due to the typical matching framework.  For major details see D'Orazio et al. (2006) or the package vignette.


}
\author{
Marcello D'Orazio

Maintainer: Marcello D'Orazio <mdo.statmatch@gmail.com>
}

\references{
D'Orazio M., Di Zio M., Scanu M. (2006) \emph{Statistical Matching, Theory and Practice}. Wiley, Chichester.
}

