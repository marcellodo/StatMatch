% --- Source file: mixed.mtc.Rd ---
\name{mixed.mtc}
\alias{mixed.mtc}
\title{Statistical Matching via Mixed Methods}

\description{
  This function implements some mixed methods to perform statistical matching between two data sources. 
}

\usage{
mixed.mtc(data.rec, data.don, match.vars, y.rec, z.don, method="ML",
           rho.yz=NULL, micro=FALSE, constr.alg="Hungarian") 
}

\arguments{

\item{data.rec}{
  	A matrix or data frame that plays the role of \emph{recipient} in the statistical matching application. This data set must contain all variables (columns) that should be used in statistical matching, i.e. the variables called by the arguments \cr \code{match.vars} and \code{y.rec}. Note that continuous variables are expected, if there are some categorical variables they are re-coded into dummies. Missing values (\code{NA}) are not allowed.
}

\item{data.don}{
   A matrix or data frame that plays the role of \emph{donor} in the statistical matching application. This data set must contain all the numeric variables (columns) that should be used in statistical matching, i.e. the variables called by the arguments \code{match.vars} and \code{z.don}. Note that continuous variables are expected, if there are some categorical variables they are re-coded into dummies. Missing values (\code{NA}) are not allowed. 
}

\item{match.vars}{
A character vector with the names of the common variables (the columns in both the data frames) to be used as matching variables (\bold{X}). 
}

\item{y.rec}{
A character vector with the name of the target variable Y that is observed only for units in \code{data.rec}. Only one continuous variable is allowed.
}

\item{z.don}{
A character vector with the name of the target variable Z that is observed only for units in \code{data.don}.  Only one continuous variable is allowed.
}

\item{method}{
A character vector that identifies the method that should be used to estimate the parameters of the regression models: Y vs. \bold{X} and Z vs. \bold{X}.  Maximum Likelihood method is used when \code{method="ML"} (default); on the contrary, when \code{method="MS"} the parameters are estimated according to approach proposed by Moriarity and Scheuren (2001 and 2003).  See Details for further information.
}

\item{rho.yz}{
A numeric value representing a guess for the correlation between the Y (\code{y.rec}) and the Z variable (\code{z.don}) that are not jointly observed.  When \code{method="MS"} then the argument \code{cor.yz} must specify the value of the correlation coefficient \eqn{\rho_{YZ}}{rho_YZ}; on the contrary, when \code{method="ML"}, it must specify the \emph{partial correlation coefficient} between Y and Z given \bold{X} (\eqn{\rho_{YZ|\bf{X}}}{rho_YZ|X}). 

By default (\code{rho.yz=NULL}).  In practice, in absence of auxiliary information concerning the correlation coefficient or the partial correlation coefficient, the statistical matching is carried out under the assumption of independence between Y and Z given \bold{X} (Conditional Independence Assumption, CIA ), i.e.  \eqn{\rho_{YZ|\bf{X}}=0}{rho_YZ|X = 0}. 
}

\item{micro}{
Logical.  When \code{micro=FALSE} (default) only the parameters' estimates are returned.  On the contrary, when \code{micro=TRUE} the function returns also \code{data.rec} filled in with the values for the variable Z.  The donors for filling in Z in \code{data.rec} are identified using a constrained distance hot deck method.  In this case, the number of units (rows) in \code{data.don} must be grater or equal to the number of units (rows) in \code{data.rec}.  See next argument and Details for further information.
}

\item{constr.alg}{
A string that has to be specified when \code{micro=TRUE}, in order to solve the transportation problem involved by the constrained distance hot deck method.  Two choices are available: \dQuote{lpSolve} and \dQuote{Hungarian}.  In the first case, \cr
\code{constr.alg="lpSolve"}, the transportation problem is solved by means of the function \code{\link[lpSolve]{lp.transport}} available in the package \pkg{lpSolve}.  When \cr \code{constr.alg="Hungarian"} (default) the transportation problem is solved using the Hungarian method implemented in the function \code{\link[clue]{solve_LSAP}} available in the package \pkg{clue} (Hornik, 2012).  Note that Hungarian algorithm is more efficient and requires less processing time.
}

}
  
\details{
This function implements some mixed methods to perform statistical matching. A mixed method consists of two steps: 

(1) adoption of a parametric model for the joint distribution of \eqn{ \left( \mathbf{X},Y,Z \right) }{(\bold{X},Y,Z)} and estimation of its parameters;

(2) derivation of a complete \dQuote{synthetic} data set (recipient data set filled in with values for the Z variable) using a nonparametric approach.


In this case, as far as (1) is concerned, it is assumed that  \eqn{ \left( \mathbf{X},Y,Z \right) }{(\bold{X},Y,Z)} follows a multivariate normal distribution.  Please note that if some of the \bold{X} are categorical, then they are recoded into dummies before starting with the estimation.  In such a case, the assumption of multivariate normal distribution may be questionable.  

The whole procedure is based on the imputation method known as \emph{predictive mean matching}.  The procedure consists of three steps: 


\bold{step 1a)} \emph{Regression step}: the two linear regression models Y vs. \bold{X} and Z vs. \bold{X} are considered and their parameters are estimated. 


\bold{step 1b)} \emph{Computation of intermediate values}.  For the units in \code{data.rec} the following intermediate values are derived:

\deqn{ \tilde{z}_{a} = \hat{\alpha}_{Z} + \hat{\beta}_{Z\bf{X}} \mathbf{x}_a + e_a }{z_a = alpha_Z + beta_ZX * x_a + e_a }

for each \eqn{a=1,\ldots,n_{A}}{a=1,...,n_A}, being \eqn{n_A}{n_A} the number of units in \code{data.rec} (rows of \code{data.rec}).  Note that, \eqn{e_a}{e_a} is a random draw from the multivariate normal distribution with zero mean and estimated residual variance  \eqn{\hat{\sigma}_{Z|\bf{X}}}{sigma_ZX}.

Similarly, for the units in \code{data.don} the following intermediate values are derived:

\deqn{ \tilde{y}_{b} = \hat{\alpha}_{Y} + \hat{\beta}_{Y\bf{X}} \mathbf{x}_b + e_b }{ y_b = alpha_Y + beta_YX * x_b + e_b  }

for each \eqn{b=1,\ldots,n_{B}}{1,...,n_B}, being \eqn{n_B}{n_B} the number of units in \code{data.don} (rows of \code{data.don}). \eqn{e_b}{e_b} is a random draw from the multivariate normal distribution with zero mean and estimated residual variance \eqn{\hat{\sigma}_{Y|\bf{X}}}{sigma_YX}.


\bold{step 2)} \emph{Matching step}.  For each observation (row) in \code{data.rec} a donor is chosen in \code{data.don} through a nearest neighbor constrained distance hot deck procedure.  The distances are computed between \eqn{\left( y_a, \tilde{z}_a \right)}{(y_a, z^_a)} and \eqn{\left( \tilde{y}_b, z_b \right)}{(y^_b, z_b)} using Mahalanobis distance.


For further details see Sections 2.5.1 and 3.6.1 in D'Orazio \emph{et al.} (2006).

In step 1a) the parameters of the regression model can be estimated by means of the Maximum Likelihood method (\code{method="ML"}) (see D'Orazio \emph{et al.}, 2006, pp. 19--23,73--75) or, using the Moriarity and Scheuren (2001 and 2003) approach (\code{method="MS"}) (see also D'Orazio \emph{et al.}, 2006, pp. 75--76).  The two estimation methods are compared in D'Orazio \emph{et al.} (2005). 

When \code{method="MS"}, if the value specified for the argument \code{rho.yz} is not compatible with the other correlation coefficients estimated from the data, then it is substituted with the closest value compatible with the other estimated coefficients.
 
When \code{micro=FALSE} only the estimation of the parameters is performed (step 1a).  Otherwise, \cr 
(\code{micro=TRUE}) the whole procedure is carried out.

} 

\value{
A list with a varying number of components depending on the values of the arguments 
\code{method} and \code{rho.yz}. 

\item{mu}{
The estimated mean vector. 
}

\item{vc}{
The estimated variance--covariance matrix. 
}

\item{cor}{
The estimated correlation matrix. 
}

\item{res.var}{
A vector with estimates of the residual variances \eqn{\sigma_{Y|Z\bf{X}}}{ sigma_Y|ZX} and \eqn{\sigma_{Z|Y\bf{X}}}{ sigma_Z|YX}. 
}

\item{start.prho.yz}{
It is the initial guess for the partial correlation coefficient \eqn{\rho_{YZ|\bf{X}}}{rho_YZ|X} passed in input via the \code{rho.yz} argument when \code{method="ML"}.
}

\item{rho.yz}{
Returned in output only when \code{method="MS"}. It is a vector with four values: the initial guess for \eqn{\rho_{YZ}}{ rho_YZ}; the lower and upper bounds for \eqn{\hat{\rho}_{YZ}}{rho_YZ} in the statistical matching framework given the correlation coefficients between Y and \bold{X} and the correlation coefficients between Z and \bold{X} estimated from the available data; and, finally, the closest admissible value used in computations instead of the initial \code{rho.yz} that resulted not coherent with the others correlation coefficients estimated from the available data.
}

\item{phi}{
When \code{method="MS"}. Estimates of the \eqn{\phi}{phi} terms introduced by Moriarity and Scheuren (2001 and 2003). 
}



\item{filled.rec}{
The \code{data.rec} filled in with the values of Z. It is returned only when \cr \code{micro=TRUE}.  
}

\item{mtc.ids}{
when \code{micro=TRUE}. This is a matrix with the same number of rows of \code{data.rec} and two columns. The first column contains the row names of the \code{data.rec} and the second column contains the row names of the corresponding donors selected from the \code{data.don}. When the input matrices do not contain row names, a numeric matrix with the indexes of the rows is provided.
}

\item{dist.rd}{
A vector with the distances between each recipient unit and the corresponding donor, returned only in case \code{micro=TRUE}.
}

\item{call}{
How the function has been called.
}

}


\references{

D'Orazio, M., Di Zio, M. and Scanu, M. (2005). \dQuote{A comparison among different estimators of regression parameters on statistically matched files through an extensive simulation study}, \emph{Contributi}, \bold{2005/10}, Istituto Nazionale di Statistica, Rome.

D'Orazio, M., Di Zio, M. and Scanu, M. (2006). \emph{Statistical Matching: Theory and Practice.} Wiley, Chichester.

Hornik K. (2012).  clue: Cluster ensembles.  R package version 0.3-45.  \url{https://CRAN.R-project.org/package=clue}.

Moriarity, C., and Scheuren, F. (2001). \dQuote{Statistical matching: a paradigm for assessing the uncertainty in the procedure}. \emph{Journal of Official Statistics}, \bold{17}, 407--422.

Moriarity, C., and Scheuren, F. (2003). \dQuote{A note on Rubin's statistical matching using file concatenation with adjusted weights and multiple imputation}, \emph{Journal of Business and Economic Statistics}, \bold{21}, 65--73.

}

\author{
 Marcello D'Orazio \email{mdo.statmatch@gmail.com} 
}

\seealso{ 
\code{\link[StatMatch]{NND.hotdeck}}, \code{\link[StatMatch]{mahalanobis.dist}}
}

\examples{

# reproduce the statistical matching framework
# starting from the iris data.frame
suppressWarnings(RNGversion("3.5.0"))
set.seed(98765)
pos <- sample(1:150, 50, replace=FALSE)
ir.A <- iris[pos,c(1,3:5)]
ir.B <- iris[-pos, 2:5]

xx <- intersect(colnames(ir.A), colnames(ir.B))
xx  # common variables

# ML estimation method under CIA ((rho_YZ|X=0));
# only parameter estimates (micro=FALSE)
# only continuous matching variables
xx.mtc <- c("Petal.Length", "Petal.Width")
mtc.1 <- mixed.mtc(data.rec=ir.A, data.don=ir.B, match.vars=xx.mtc,
                    y.rec="Sepal.Length", z.don="Sepal.Width")

# estimated correlation matrix
mtc.1$cor 

# ML estimation method under CIA ((rho_YZ|X=0));
# only parameter estimates (micro=FALSE)
# categorical variable 'Species' used as matching variable

xx.mtc <- xx
mtc.2 <- mixed.mtc(data.rec=ir.A, data.don=ir.B, match.vars=xx.mtc,
                    y.rec="Sepal.Length", z.don="Sepal.Width")

# estimated correlation matrix
mtc.2$cor 


# ML estimation method with partial correlation coefficient
# set equal to 0.5 (rho_YZ|X=0.5)
# only parameter estimates (micro=FALSE)

mtc.3 <- mixed.mtc(data.rec=ir.A, data.don=ir.B, match.vars=xx.mtc,
                    y.rec="Sepal.Length", z.don="Sepal.Width",
                    rho.yz=0.5)

# estimated correlation matrix
mtc.3$cor 

# ML estimation method with partial correlation coefficient
# set equal to 0.5 (rho_YZ|X=0.5)
# with imputation step (micro=TRUE)

mtc.4 <- mixed.mtc(data.rec=ir.A, data.don=ir.B, match.vars=xx.mtc,
                    y.rec="Sepal.Length", z.don="Sepal.Width",
                    rho.yz=0.5, micro=TRUE, constr.alg="Hungarian")

# first rows of data.rec filled in with z
head(mtc.4$filled.rec)

#
# Moriarity and Scheuren estimation method under CIA;
# only with parameter estimates (micro=FALSE)
mtc.5 <- mixed.mtc(data.rec=ir.A, data.don=ir.B, match.vars=xx.mtc,
                    y.rec="Sepal.Length", z.don="Sepal.Width",
                    method="MS")

# the starting value of rho.yz and the value used
# in computations
mtc.5$rho.yz

# estimated correlation matrix
mtc.5$cor 

# Moriarity and Scheuren estimation method
# with correlation coefficient set equal to -0.15 (rho_YZ=-0.15)
# with imputation step (micro=TRUE)

mtc.6 <- mixed.mtc(data.rec=ir.A, data.don=ir.B, match.vars=xx.mtc,
                    y.rec="Sepal.Length", z.don="Sepal.Width",
                    method="MS", rho.yz=-0.15, 
                    micro=TRUE, constr.alg="lpSolve")

# the starting value of rho.yz and the value used
# in computations
mtc.6$rho.yz

# estimated correlation matrix
mtc.6$cor

# first rows of data.rec filled in with z imputed values
head(mtc.6$filled.rec)

}

\keyword{nonparametric}
\keyword{regression}