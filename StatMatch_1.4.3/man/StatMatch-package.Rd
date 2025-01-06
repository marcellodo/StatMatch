\name{StatMatch-package}
\alias{StatMatch-package}
\alias{StatMatch}
\docType{package}
\title{
Statistical Matching or Data Fusion
}
\description{
Functions to perform statistical matching (aka data fusion), i.e. the integration of two data sources.  Some functions can also be used to impute missing values in data sets through hotdeck imputation methods.

}
\details{

\emph{Statistical matching} (also known as \emph{data fusion}) aims to integrate two data sources (typically samples) that refer to the same target population and share a number of variables.  The ultimate aim is to study the relationship between variables that are not jointly observed in a single data source.  The integration can be performed at \emph{micro} level (a \emph{synthetic} or \emph{fused} file is the output) or \emph{macro} level (the aim is to estimate correlation coefficients, regression coefficients, contingency tables, etc.).  
Non-parametric hotdeck imputation methods (random, rank and nearest neighbour) can be used to derive the synthetic data set.  Alternatively, a mixed (parametric-nonparametric) procedure based on predictive mean matching can be used.  Methods are also available for performing statistical matching when dealing with data from complex sample surveys.  Finally, some functions can be used to explore the uncertainty due to the typical matching framework and can help to decide whether or not to proceed with statistical matching at the micro or macro level.  For major details see D'Orazio et al. (2006), D'Orazio (2015), and material available at \url{https://github.com/marcellodo/StatMatch/tree/master/Tutorials_Vignette_OtherDocs}.


}
\author{
Marcello D'Orazio

Maintainer: Marcello D'Orazio <mdo.statmatch@gmail.com>
}

\references{
D'Orazio M. (2015) Integration and imputation of survey data in R: the StatMatch package, \emph{Romanian Statistical Review}, 2/2015, pp. 57--68

D'Orazio M., Di Zio M., Scanu M. (2006) \emph{Statistical Matching, Theory and Practice}. Wiley, Chichester.
}

