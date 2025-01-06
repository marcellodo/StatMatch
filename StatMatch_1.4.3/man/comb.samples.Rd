% --- Source file:  ---
\name{comb.samples}
\alias{comb.samples}

\title{Statistical Matching of data from complex sample surveys}

\description{
This function allows you to cross-tabulate two categorical variables, Y and Z, observed separately in two independent surveys (Y is collected in survey A and Z is collected in survey B) relating to the same target population.  The two surveys share a number of common variables \bold{X}.  If a third survey C is available on the same population and collects both Y and Z, these data are used as a source of additional information.

The statistical adjustment is done by calibrating the survey weights as suggested in Renssen (1998).

It is also possible to use the function to derive estimates that a unit falls into one of the categories of the target variable (the estimates are based on Linear probability models and are obtained as a by-product of the Renssen method).
}

\usage{
comb.samples(svy.A, svy.B, svy.C=NULL, y.lab, z.lab, form.x, 
              estimation=NULL, micro=FALSE, ...)
}

\arguments{

\item{svy.A}{
A \code{svydesign} \R object that stores the data collected in the survey A and all the information concerning the corresponding sampling design.  This object can be created by using the function \code{\link[survey]{svydesign}} in the package \pkg{survey}.  All the variables specified in \code{form.x} and by \code{y.lab} must be available in \code{svy.A}.
}

\item{svy.B}{
A \code{svydesign} \R object that stores the data collected in the survey B and all the information concerning the corresponding sampling design.  This object can be created by using the function \code{\link[survey]{svydesign}} in the package \pkg{survey}.  All the variables specified in \code{form.x} and by \code{z.lab} must be available in \code{svy.B}.
}

\item{svy.C}{
A \code{svydesign} \R object that stores the data collected in the the survey C and all the information concerning the corresponding sampling design.  This object can be created by using the function \code{\link[survey]{svydesign}} in the package \pkg{survey}.

When \code{svy.C=NULL} (default), i.e. no auxiliary information is available, the function returns an estimate of the contingency table of Y vs. Z under the Conditional Independence assumption (CIA) (see Details for major information).

When \code{svy.C} is available, if \code{estimation="incomplete"} then it must contain at least \code{y.lab} and \code{z.lab} variables.  On the contrary, when 
\code{estimation="synthetic"} all the variables specified in \code{form.x}, \code{y.lab} and \code{z.lab} must be available in \code{svy.C}.

}

\item{y.lab}{
A string providing the name of the Y variable, available in survey A and in survey C (if available).  The Y variable can be a categorical variable (\code{factor} in \R) or a continuous one (in this latter case \code{z.lab} should be categorical).
}

\item{z.lab}{
A string providing the name of the Z variable available in survey B and in survey C (if available).  The Z variable can be a categorical variable (\code{factor} in \R) or a continuous one (in this latter case \code{y.lab} should be categorical).
}


\item{form.x}{
A \R formula specifying which matching variables (subset of the \bold{X} variables) are collected in all surveys and how they are to be considered when combining samples.  For example, \code{form.x=~x1+x2} means that variables x1 and x2 are to be considered marginally, without taking into account their cross-tabulation; only their marginal distribution is considered.  To skip the intercept, the formula must be written as \code{form.x=~x1+x2-1}.

When dealing with categorical variables, \code{form.x=~x1:x2-1} means that the joint distribution of the two variables (table of x1 vs. x2) has to be considered.

To better understand the use of \code{form.x}, see \code{\link[stats]{model.matrix}} (see also \cr
\code{\link[stats]{formula}}).

Due to the weight calibration features, it is preferable to work with categorical \bold{X} variables.  In some cases, the procedure can be successful when a single continuous variable is considered together with one or more categorical variables, but often it may be necessary to categorise the continuous variable (see Details).
}

\item{estimation}{
A string identifying the method to be used to estimate the table of Y vs. Z when data from survey C are available.  As suggested in Renssen (1998), there are two alternative methods: (i) incomplete two-way stratification (\code{estimation="incomplete"} or \code{estimation="ITWS"}, the default) and (ii) synthetic two-way stratification (\code{estimation="synthetic"} or \cr \code{estimation="STWS"}).  In the first case (\code{estimation="incomplete"}), only Y and Z variables need to be available in \code{svy.C}.  On the other hand, in the case of \cr \code{estimation="synthetic"}, the C survey must contain all the \bold{X} variables specified by \code{form.x}, as well as the Y and Z variables.  See details for more information.
}

\item{micro}{
Logical, when \code{TRUE} predictions of Z in A and of Y in B are provided.  In particular when Y and Z are both categorical variables it is provided the estimated probability that a unit falls in each of the  categories of the given variable.  These probabilities are estimated as a by-product of the whole procedure by considering Linear Probability Models, as suggested in Renssen (1998) (see Details)
}

\item{...}{

Further arguments that may be necessary for calibration.  In particular, the argument \code{calfun} allows to specify the calibration function: \cr
(i) \code{calfun="linear"} for linear calibration (default); \cr
(ii) \code{calfun="raking"} to rake the survey weights; and \cr
(iii) \code{calfun="logit"} for logit calibration.  See \code{\link[survey]{calibrate}} for major details.

Note that when \code{calfun="linear"} calibration may return negative weights.  Generally speaking, in sample surveys weights are expected to be greater than or equal to 1, i.e. \code{bounds=c(1, Inf)}.

The number of iterations used in calibration can be modified by using the argument \code{maxit} (by default \code{maxit=50}).

See \code{\link[survey]{calibrate}} for further details.

}

}

\details{
TThis function estimates the Y vs. Z contingency table by performing a series of calibrations of the survey weights.  In practice, the estimation is carried out on data from survey C, using all the information from surveys A and B. If survey C is not available, the table of Y vs. Z is estimated under the assumption of conditional independence (CIA), i.e. \eqn{p(Y,Z)=p(Y|\bold{X}) \times p(Z|\bold{X}) \times p(\bold{X})}{p(Y,X)=p(Y|X)*p(Z|X)*p(X)}.

When data from survey C are available (Renssen, 1998), the table of Y vs. Z can be estimated by: Incomplete Two-Way Stratification (ITWS) or Synthetic Two-Way Stratification (STWS).  In the first case (ITWS) the weights of the units in survey C are calibrated so that the new weights allow to reproduce the marginal distributions of Y estimated on survey A, and that of Z estimated on survey B.  Note that the distribution of the \bold{X} variables in survey A and in survey B, must be harmonized before performing ITWS (see \code{\link[StatMatch]{harmonize.x}}).  

The Synthetic Two-Way Stratification allows to estimate the table of Y vs. Z by considering also the \bold{X} variables observed in C.  This method consists in correcting the table of Y vs. Z estimated under the CIA according to the relationship between Y and Z observed in survey C (for major details see Renssen, 1998.

If the argument \code{micro} is set to \code{TRUE} the function will also return \code{Z.A} and \code{Y.B}.  The first data.frame has the same rows as \code{svy.A} and the number of columns is equal to the number of categories of the Z variable specified by \code{z.lab}.  Each row gives the estimated probabilities of taking a value in each category.  The same happens for \code{Y.B}, which gives the estimated probabilities of assuming a category of \code{y.lab} for each unit in B. The estimated probabilities are obtained by applying the linear probability models (for more details see Renssen, 1998).  Unfortunately, such models may give estimated probabilities less than 0 or greater than 1. Great care should be taken when using such predictions for practical purposes.
}


\value{
A \R list with the results of the calibration procedure according to the input arguments. 

\item{yz.CIA}{
The table of Y (\code{y.lab}) vs. Z (\code{z.lab}) estimated under the Conditional Independence Assumption (CIA).
}

\item{cal.C}{
The survey object \code{svy.C} after the calibration.  Only when \code{svy.C} is provided.
}

\item{yz.est}{
The table of Y (\code{y.lab}) vs. Z (\code{z.lab}) estimated under the method specified via \code{estimation} argument.  Only when \code{svy.C} is provided.
}

\item{Z.A}{
Only when \code{micro=TRUE}.  It is a data frame with the same rows as in \code{svy.A} and the number of columns is equal to the number of categories of the variable \code{z.lab}. Each row provides the estimated probabilities for a unit being in the various categories of \code{z.lab}.
}

\item{Y.B}{
Only when \code{micro=TRUE}.  It is a data frame with the same rows as in \code{svy.B} and the number of columns is equal to the number of categories of \code{y.lab}. Each row provides the estimated probabilities for a unit being in the various categories of \code{y.lab}.
}

\item{call}{
Stores the call to this function with all the values specified for the various arguments (\code{call=match.call()}).
}

}

\references{

D'Orazio, M., Di Zio, M. and Scanu, M. (2006). \emph{Statistical Matching: Theory and Practice}. Wiley, Chichester.

Renssen, R.H. (1998) \dQuote{Use of Statistical Matching Techniques in Calibration Estimation}. \emph{Survey Methodology}, \bold{24}, pp. 171--183.

}


\author{
 Marcello D'Orazio \email{mdo.statmatch@gmail.com} 
}

\seealso{ 
\code{\link[survey]{calibrate}}, \code{\link[survey]{svydesign}}, \code{\link[StatMatch]{harmonize.x}}
}

\examples{

data(quine, package="MASS") #loads quine from MASS
str(quine)
quine$c.Days <- cut(quine$Days, c(-1, seq(0,20,10),100))
table(quine$c.Days)


# split quine in two subsets

suppressWarnings(RNGversion("3.5.0"))
set.seed(124)
lab.A <- sample(nrow(quine), 70, replace=TRUE)
quine.A <- quine[lab.A, c("Eth","Sex","Age","Lrn")]
quine.B <- quine[-lab.A, c("Eth","Sex","Age","c.Days")]

# create svydesign objects
require(survey)
quine.A$f <- 70/nrow(quine) # sampling fraction
quine.B$f <- (nrow(quine)-70)/nrow(quine)
svy.qA <- svydesign(~1, fpc=~f, data=quine.A)
svy.qB <- svydesign(~1, fpc=~f, data=quine.B)

# Harmonizazion wrt the joint distribution
# of ('Sex' x 'Age' x 'Eth')

# vector of population total known
# estimated from the full data set
# note the formula!
tot.m <- colSums(model.matrix(~Eth:Sex:Age-1, data=quine))
tot.m

out.hz <- harmonize.x(svy.A=svy.qA, svy.B=svy.qB, x.tot=tot.m,
            form.x=~Eth:Sex:Age-1, cal.method="linear")
            
# estimation of 'Lrn' vs. 'c.Days' under the CIA

svy.qA.h <- out.hz$cal.A
svy.qB.h <- out.hz$cal.B

out.1 <- comb.samples(svy.A=svy.qA.h, svy.B=svy.qB.h,
            svy.C=NULL, y.lab="Lrn", z.lab="c.Days",
            form.x=~Eth:Sex:Age-1)

out.1$yz.CIA
addmargins(out.1$yz.CIA)

#
# incomplete two-way stratification

# select a sample C from quine
# and define a survey object

suppressWarnings(RNGversion("3.5.0"))
set.seed(4321)
lab.C <- sample(nrow(quine), 50, replace=TRUE)
quine.C <- quine[lab.C, c("Lrn","c.Days")]
quine.C$f <- 50/nrow(quine) # sampling fraction
svy.qC <- svydesign(~1, fpc=~f, data=quine.C)

# call comb.samples
out.2 <- comb.samples(svy.A=svy.qA.h, svy.B=svy.qB.h,
            svy.C=svy.qC, y.lab="Lrn", z.lab="c.Days",
            form.x=~Eth:Sex:Age-1, estimation="incomplete",
            calfun="linear", maxit=100)

summary(weights(out.2$cal.C))
out.2$yz.est # estimated table of 'Lrn' vs. 'c.Days'
# difference wrt the table 'Lrn' vs. 'c.Days' under CIA
addmargins(out.2$yz.est)-addmargins(out.2$yz.CIA)

# synthetic two-way stratification
# only macro estimation

quine.C <- quine[lab.C, ]
quine.C$f <- 50/nrow(quine) # sampling fraction
svy.qC <- svydesign(~1, fpc=~f, data=quine.C)

out.3 <- comb.samples(svy.A=svy.qA.h, svy.B=svy.qB.h,
            svy.C=svy.qC, y.lab="Lrn", z.lab="c.Days",
            form.x=~Eth:Sex:Age-1, estimation="synthetic",
            calfun="linear",bounds=c(.5,Inf), maxit=100)

summary(weights(out.3$cal.C))

out.3$yz.est # estimated table of 'Lrn' vs. 'c.Days'
# difference wrt the table of 'Lrn' vs. 'c.Days' under CIA
addmargins(out.3$yz.est)-addmargins(out.3$yz.CIA)
# diff wrt the table of 'Lrn' vs. 'c.Days' under incomplete 2ws
addmargins(out.3$yz.est)-addmargins(out.2$yz.CIA)

# synthetic two-way stratification
# with micro predictions

out.4 <- comb.samples(svy.A=svy.qA.h, svy.B=svy.qB.h,
            svy.C=svy.qC, y.lab="Lrn", z.lab="c.Days",
            form.x=~Eth:Sex:Age-1, estimation="synthetic",
            micro=TRUE, calfun="linear",bounds=c(.5,Inf), 
            maxit=100)
            
head(out.4$Z.A)
head(out.4$Y.B)



}

\keyword{survey}
