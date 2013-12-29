% --- Source file:  ---
\name{harmonize.x}
\alias{harmonize.x}
\title{Harmonizes the marginal (joint) distribution of a set of variables observed independently in two sample surveys referred to the same target population}

\description{
This function permits to harmonize the marginal or the joint distribution of a set of variables observed independently in two sample surveys carried out on the same target population.  This harmonization is carried out by using the calibration of the survey weights of the sample units in both the surveys according to the procedure suggested by Renssen (1998).
}

\usage{
harmonize.x(svy.A, svy.B, form.x, x.tot=NULL, 
                      cal.method="linear", ...)
}

\arguments{

\item{svy.A}{
A \code{svydesign} \R object that stores the data collected in the the survey A and all the information concerning the sampling design.  This object can be created by using the function \code{\link[survey]{svydesign}} in the package \pkg{survey}.
}

\item{svy.B}{
A \code{svydesign} \R object that stores the data collected in the the survey B and all the information concerning the sampling design.  This object can be created by using the function \code{\link[survey]{svydesign}} in the package \pkg{survey}.}

\item{form.x}{
A \R formula specifying which of the variables, common to both the surveys, have to be considered, and how have to be considered.  For instance \cr
\code{form.x=~x1+x2} means that the marginal distribution of the variables x1 and x2 have to be harmonized and there is also an \sQuote{Intercept}.  In order to skip the intercept the formula has to be written in the following manner \cr \code{form.x=~x1+x2-1}.

When dealing with categorical variables, if the formula is written in the following manner \code{form.x=~x1:x2-1} it means that the harmonization has to be carried out by considering the joint distribution of the two variables (x1 vs. x2).  To better understand how \code{form.x} works see \code{\link[stats]{model.matrix}} (see also \code{\link[stats]{formula}}).

Due to weights calibration features, it is preferable to work with categorical \bold{X} variables.  In some cases, the procedure may be successful when a single continuous variable is considered jointly with one or more categorical variables.  When dealing with several continuous variable it may be preferable to categorize them.
}

\item{x.tot}{
A vector or table with known population totals for the \bold{X} variables.  A vector is required when \code{cal.method="linear"} or \code{cal.method="raking"}.  The names and the length of the vector depends on the way it is specified the argument \code{form.x} (see \code{\link[stats]{model.matrix}}).  A contingency table is required when \cr \code{cal.method="poststratify"} (for details see \code{\link[survey]{postStratify}}).

When \code{x.tot} is not provided (i.e. \code{x.tot=NULL}) then the vector of totals is estimated as a weighted average of the totals estimated on the two surveys.  The weight assigned to the totals estimated from A is \eqn{\lambda = n_A/(n_A+n_B)}{lambda= n_A/(n_A+n_B)};  \eqn{1-\lambda}{1-lambda} is the weight assigned to \bold{X} totals estimated from survey B (\eqn{n_A}{n_A} and \eqn{n_B}{n_B} are the number of units in A and B respectively).
}

\item{cal.method}{
A string that specifies how the calibration of the weights in \code{svy.A} and \code{svy.B} has to be carried out.  By default \code{cal.method="linear"} linear calibration is performed.  In particular, the calibration is carried out by mean of the function \code{\link[survey]{calibrate}} in the package \pkg{survey}.

Alternatively, it is possible to rake the origin survey weights by specifying \code{cal.method="raking"}.  Finally, it is possible to perform a simple post-stratification by setting \code{cal.method="poststratify"}.  Note that in this case the weights adjustments are carried out by considering the function \cr
\code{\link[survey]{postStratify}} in the package \pkg{survey}.
}
\item{...}{

Further arguments that may be necessary for calibration or post-stratification.  In particular, when \code{cal.method="linear"} there is the chance of having negative weights.  This drawback can be avoided by requiring that the final weights lie in a given interval specified via \code{bounds} argument (an argument of \code{\link[survey]{calibrate}}).  In sample surveys  one expects that weights are greater than or equal to 1, i.e. \code{bounds=c(1, Inf)}.

The number of iterations used in calibration can be modified too by using the argument \code{maxit} (by default \code{maxit=50}).

See \code{\link[survey]{calibrate}} for further details.

}

}

\details{
This function harmonizes the totals of the \bold{X} variables, observed in both survey A and survey B, to be equal to given known totals specified via \code{x.tot}.  When these totals are not known (\code{x.tot=NULL})  they are estimated by combining the estimates derived from the two separate surveys.  The harmonization is carried out according to a procedure suggested by Renssen (1998) based on calibration of survey weights (for major details on calibration see Sarndal and Lundstrom, 2005).  The procedure is particularly suited to deal with categorical \bold{X} variables, in this case it permits to harmonize the joint or the marginal distribution of the categorical variables being considered. Note that an incomplete crossing of the \bold{X} variables can be considered: i.e. harmonisation wrt to the joint distribution of \eqn{X_1 \times X_2}{X_1 x X_2} and the marginal distribution of \eqn{X_3}{x_3}).

The calibration procedure may not produce the final result due to convergence problems.  In this case an error message appears.  In order to reach convergence it may be necessary to launch the procedure with less constraints (i.e a reduced number of population totals) by joining adjacent categories or by discarding some variables.  

In some limited cases it could be possible to consider both categorical and continuous variables.  In this situation it may happen that calibration is not successful.  In order to reach convergence it may be necessary to categorize the continuous \bold{X} variables.

Post-stratification is a special case of calibration; all the weights of the units in a given post-stratum are modified so the reproduce the known total for that post-stratum.  Post-stratification avoids problems of convergence but, on the other hand it may produce final weights with a higher variability than those derived from the calibration.

}


\value{
A \R with list the results of calibration procedures carried out on survey A and survey B, respectively.  In particular the following components will be provided:

\item{cal.A}{
The survey object \code{svy.A} after the calibration; in particular, the weights now are calibrated with respect to the totals of the \bold{X} variables.
}

\item{cal.B}{
The survey object \code{svy.B} after the calibration; in particular, the weights now are calibrated with respect to the totals of the \bold{X} variables.
}

\item{weights.A}{
The new calibrated weights associated to the the units in \code{svy.A}.
}

\item{weights.B}{
The new calibrated weights associated to the the units in \code{svy.B}.
}

\item{call}{
Stores the call to this function  with all the values specified for the various arguments (\code{call=match.call()}).
}

}

\references{

D'Orazio, M., Di Zio, M. and Scanu, M. (2006). \emph{Statistical Matching: Theory and Practice}. Wiley, Chichester.

Renssen, R.H. (1998) \dQuote{Use of Statistical Matching Techniques in Calibration Estimation}. \emph{Survey Methodology}, N. \bold{24}, pp. 171--183.

Sarndal, C.E. and Lundstrom, S. (2005) \emph{Estimation in Surveys with Nonresponse}. Wiley, Chichester.

}


\author{
 Marcello D'Orazio \email{madorazi@istat.it} 
}

\seealso{ 
\code{\link[StatMatch]{comb.samples}}, \code{\link[survey]{calibrate}}, \code{\link[survey]{svydesign}}, \code{\link[survey]{postStratify}}, 

}

\examples{

data(quine, package="MASS") #loads quine from MASS
str(quine)

# split quine in two subsets
set.seed(7654)
lab.A <- sample(nrow(quine), 70, replace=TRUE)
quine.A <- quine[lab.A, c("Eth","Sex","Age","Lrn")]
quine.B <- quine[-lab.A, c("Eth","Sex","Age","Days")]

# create svydesign objects
require(survey)
quine.A$f <- 70/nrow(quine) # sampling fraction
quine.B$f <- (nrow(quine)-70)/nrow(quine)
svy.qA <- svydesign(~1, fpc=~f, data=quine.A)
svy.qB <- svydesign(~1, fpc=~f, data=quine.B)

#------------------------------------------------------
# example (1)
# Harmonizazion of the distr. of Sex vs. Age
# usign poststratification

# (1.a) known population totals
# the population toatal are computed on the full data frame
tot.sex.age <- xtabs(~Sex+Age, data=quine)
tot.sex.age

out.hz <- harmonize.x(svy.A=svy.qA, svy.B=svy.qB, form.x=~Sex+Age,
          x.tot=tot.sex.age, cal.method="poststratify")

tot.A <- xtabs(out.hz$weights.A~Sex+Age, data=quine.A)
tot.B <- xtabs(out.hz$weights.B~Sex+Age, data=quine.B)

tot.sex.age-tot.A
tot.sex.age-tot.B

# (1.b) unknown population totals (x.tot=NULL)
# the population total is estimated by combining totals from the
 # two surveys

out.hz <- harmonize.x(svy.A=svy.qA, svy.B=svy.qB, form.x=~Sex+Age,
          x.tot=NULL, cal.method="poststratify")

tot.A <- xtabs(out.hz$weights.A~Sex+Age, data=quine.A)
tot.B <- xtabs(out.hz$weights.B~Sex+Age, data=quine.B)

tot.A
tot.A-tot.B

#-----------------------------------------------------
# example (2)
# Harmonizazion wrt the maginal distribution
# of 'Eth', 'Sex' and 'Age'
# using linear calibration

# (2.a) vector of population total known
# estimated from the full data set
# note the formula! only marginal distribution of the
# variables are considered
tot.m <- colSums(model.matrix(~Eth+Sex+Age-1, data=quine))
tot.m

out.hz <- harmonize.x(svy.A=svy.qA, svy.B=svy.qB, x.tot=tot.m,
            form.x=~Eth+Sex+Age-1, cal.method="linear")

summary(out.hz$weights.A) #check for negative weights
summary(out.hz$weights.B) #check for negative weights

tot.m
svytable(formula=~Eth, design=out.hz$cal.A)
svytable(formula=~Eth, design=out.hz$cal.B)

svytable(formula=~Sex, design=out.hz$cal.A)
svytable(formula=~Sex, design=out.hz$cal.B)

# Note: margins are equal but joint distributions are not!
svytable(formula=~Sex+Age, design=out.hz$cal.A)
svytable(formula=~Sex+Age, design=out.hz$cal.B)

# (2.b) vector of population total unknown
out.hz <- harmonize.x(svy.A=svy.qA, svy.B=svy.qB, x.tot=NULL,
            form.x=~Eth+Sex+Age-1, cal.method="linear")
svytable(formula=~Eth, design=out.hz$cal.A)
svytable(formula=~Eth, design=out.hz$cal.B)

svytable(formula=~Sex, design=out.hz$cal.A)
svytable(formula=~Sex, design=out.hz$cal.B)

#-----------------------------------------------------
# example (3)
# Harmonizazion wrt the joint distribution of 'Sex' vs. 'Age'
# and the marginal distribution of 'Eth'
# using raking

# vector of population total known
# estimated from the full data set
# note the formula!
tot.m <- colSums(model.matrix(~Eth+(Sex:Age-1)-1, data=quine))
tot.m

out.hz <- harmonize.x(svy.A=svy.qA, svy.B=svy.qB, x.tot=tot.m,
            form.x=~Eth+(Sex:Age)-1, cal.method="raking")

summary(out.hz$weights.A) #check for negative weights
summary(out.hz$weights.B) #check for negative weights

tot.m
svytable(formula=~Eth, design=out.hz$cal.A)
svytable(formula=~Eth, design=out.hz$cal.B)

svytable(formula=~Sex+Age, design=out.hz$cal.A)
svytable(formula=~Sex+Age, design=out.hz$cal.B)

#-----------------------------------------------------
# example (4)
# Harmonizazion wrt the joint distribution
# of ('Sex' x 'Age' x 'Eth')

# vector of population total known
# estimated from the full data set
# note the formula!
tot.m <- colSums(model.matrix(~Eth:Sex:Age-1, data=quine))
tot.m

out.hz <- harmonize.x(svy.A=svy.qA, svy.B=svy.qB, x.tot=tot.m,
            form.x=~Eth:Sex:Age-1, cal.method="linear")
tot.m
svytable(formula=~Eth+Sex+Age, design=out.hz$cal.A)
svytable(formula=~Eth+Sex+Age, design=out.hz$cal.B)

}

\keyword{survey}
