% --- Source file: gower_dist.Rd ---
\name{gower.dist}
\alias{gower.dist}
\title{Computes the Gower's Distance}

\description{
  This function computes the Gower's distance (dissimilarity) among units in a dataset or among observations in two distinct datasets.
}

\usage{
gower.dist(data.x, data.y=data.x, rngs=NULL, KR.corr=TRUE)
}

\arguments{

\item{data.x}{
  	A matrix or a data frame containing variables that should be used in the computation of the distance.
    
    Columns of mode \code{numeric} will be considered as interval scaled variables; columns of mode \code{character} or class \code{factor} will be considered as categorical nominal variables; columns of class \code{ordered} will be considered as categorical ordinal variables and, columns of mode \code{logical} will be considered as binary asymmetric variables (see Details for further information). 

Missing values (\code{NA}) are allowed.

If only \code{data.x} is supplied, the dissimilarities between rows of \code{data.x} will be computed. 
}

\item{data.y}{
   A numeric matrix or data frame with the same variables, of the same type, as those in \code{data.x}. Dissimilarities between rows of \code{data.x} and rows of \code{data.y} will be computed. If not provided, by default it is assumed equal to \code{data.x} and only dissimilarities between rows of \code{data.x} will be computed. 
}

\item{rngs}{
A vector with the ranges to scale the variables. Its length must be equal to number of variables in \code{data.x}. In correspondence of nonnumeric variables, just put 1 or \code{NA}. When \code{rngs=NULL} (default) the range of a numeric variable is estimated by jointly considering the values for the variable in \code{data.x} and those in \code{data.y}. Therefore, assuming \code{rngs=NULL}, if a variable \code{"X1"} is considered:
\preformatted{
rngs["X1"] <- max(data.x[,"X1"], data.y[,"X1"]) - 
               min(data.x[,"X1"], data.y[,"X1"])}.
}

\item{KR.corr}{
When \code{TRUE} (default) the extension of the Gower's dissimilarity measure proposed by Kaufman and Rousseeuw (1990) is used. Otherwise, when \cr
\code{KR.corr=FALSE}, the Gower's (1971) formula is considered.
}
}
  
\details{
This function computes distances among records when variables of different type (categorical and continuous) have been observed. In order to handle different types of variables, the Gower's  dissimilarity coefficient (Gower, 1971) is used. By default (\code{KR.corr=TRUE}) the Kaufman and Rousseeuw (1990) extension of the Gower's dissimilarity coefficient is used. 

The final dissimilarity between the \emph{i}th and \emph{j}th unit is obtained as a weighted sum of dissimilarities for each variable:
\deqn{d(i,j) = \frac{\sum_k{\delta_{ijk} d_{ijk}}}{\sum_k{\delta_{ijk}}}}{d(i,j) = sum_k(delta_ijk * d_ijk ) / sum_k( delta_ijk )} 

In particular, \eqn{d_{ijk}}{d_ijk} represents the distance between the \emph{i}th and \emph{j}th unit computed considering the \emph{k}th variable. It depends on the nature of the variable:
\itemize{

		\item \code{logical} columns are considered as asymmetric binary variables, for such case \eqn{d_{ijk}=0}{d_ijk = 0} if \eqn{x_{ik} = x_{jk} = \code{TRUE}}{x_ik = x_jk = TRUE}, 1 otherwise;
		
		\item \code{factor} or \code{character} columns are considered as categorical nominal variables and \eqn{d_{ijk}=0}{d_ijk = 0} if \eqn{x_{ik}=x_{jk}}{x_ik = x_jk}, 1 otherwise;
		
		\item \code{numeric} columns are considered as interval-scaled variables and 
    \deqn{d_{ijk}=\frac{\left|x_{ik}-x_{jk}\right|}{R_k}}{d_ijk = abs(x_ik - x_jk) / R_k}
    being  \eqn{R_k}{R_k} the range of the \emph{k}th variable. The range is the one supplied with the argument \code{rngs}  (\code{rngs[k]}) or the one computed on available data (when \code{rngs=NULL});
		
		\item \code{ordered} columns are considered as categorical ordinal variables and the values are substituted with the corresponding position index, \eqn{r_{ik}}{r_ik} in the factor levels. When \code{KR.corr=FALSE} these position indexes (that are different from the output of the R function \code{\link[base]{rank}}) are transformed in the following manner  
    \deqn{z_{ik}=\frac{(r_{ik}-1)}{max\left(r_{ik}\right) - 1}}{z_ik = (r_ik - 1)/(max(r_ik) - 1)}
    These new values, \eqn{z_{ik}}{z_ik}, are treated as observations of an interval scaled variable.

}

As far as the weight \eqn{\delta_{ijk}}{delta_ijk} is concerned:
\itemize{

	\item \eqn{\delta_{ijk}=0}{delta_ijk = 0} if \eqn{x_{ik} = \code{NA}}{x_ik = \code{NA}} or \eqn{x_{jk} = \code{NA}}{x_jk = NA};
	
	\item \eqn{\delta_{ijk}=0}{delta_ijk = 0} if the variable is asymmetric binary  and  \eqn{x_{ik}=x_{jk}=0}{x_ik = x_jk = 0} or  \eqn{x_{ik} = x_{jk} = \code{FALSE}}{x_ik = x_jk = FALSE};
	\item 
  \eqn{\delta_{ijk}=1}{delta_ijk = 1} in all the other cases.
}

In practice, \code{NAs} and couple of cases with \eqn{x_{ik}=x_{jk}=\code{FALSE}}{x_ik = x_jk = \code{FALSE}} do not contribute to distance computation.
}

\value{
A \code{matrix} object with distances among rows of \code{data.x} and those of \code{data.y}.
}

\references{

Gower, J. C. (1971), \dQuote{A general coefficient of similarity and some of its properties}. \emph{Biometrics}, \bold{27}, 623--637.

Kaufman, L. and Rousseeuw, P.J. (1990), \emph{Finding Groups in Data: An Introduction to Cluster Analysis.} Wiley, New York.
}

\author{
 Marcello D'Orazio \email{madorazi@istat.it} 
}

\seealso{ 
\code{\link[cluster]{daisy}}, 
\code{\link[proxy]{dist}}
}

\examples{

x1 <- as.logical(rbinom(10,1,0.5)) 
x2 <- sample(letters, 10, replace=TRUE)
x3 <- rnorm(10)
x4 <- ordered(cut(x3, -4:4, include.lowest=TRUE))
xx <- data.frame(x1, x2, x3, x4, stringsAsFactors = FALSE)

# matrix of distances among observations in xx
gower.dist(xx)

# matrix of distances among first obs. in xx
# and the remaining ones
gower.dist(data.x=xx[1:3,], data.y=xx[4:10,])

}

\keyword{multivariate}
\keyword{cluster}