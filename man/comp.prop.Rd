\name{comp.prop}
\alias{comp.prop}
\title{Compares two distributions of categorical variables}

\description{
  This function compares two distributions of the same categorical variable(s).
}

\usage{
comp.prop(p1, p2, n1, n2=NULL, ref=FALSE)
}

\arguments{

\item{p1}{
 A vector or an array containing relative or absolute frequencies for one or more categorical variables.  Usually it is the output of the function \code{\link[stats]{xtabs}} or \code{\link[base]{table}}. 
}

\item{p2}{ 
A vector or an array containing relative or absolute frequencies for one or more categorical variables.  Usually it is the output of the function \code{\link[stats]{xtabs}} or \code{\link[base]{table}}.  If \code{ref=FALSE} then \code{p2} is a further estimate of the distribution of the categorical variable(s) being considered.  On the contrary (\code{ref=FALSE}) it is the 'reference' distribution (the distribution considered true).  
}

\item{n1}{
 The sample size of the sample of the sample on which \code{p1} has been estimated.
}

\item{n2}{
 The sample size of the sample of the sample on which \code{p2} has been estimated, required just when \code{ref=FALSE} (\code{p2} is estimated on another sample and is not the reference distribution).
}

\item{ref}{
	Logical.  Denotes whether \code{p2} is the reference distribution (true distribution when \code{ref=TRUE}) or just another estimate (\code{ref=FALSE}) of the distribution derived from another sample with sample size \code{n2}.
}

}


\details{
This function computes some similarity or dissimilarity measures among marginal (joint) distribution of categorical variables(s).
The following measures are considered:


\emph{Dissimilarity index} or \emph{total variation distance}: 

\deqn{\Delta_{12} = \frac{1}{2} \sum_{j=1}^J \left| p_{1,j} - p_{2,j}  \right|}{D = (1/2) * sum_j |p_1,j - p_2,j| }

where \eqn{p_{s,j}}{p_s,j} are the relative frequencies (\eqn{0 \leq p_{s,j} \leq 1}{0 <= p_s,j <= 1}). The dissimilarity index ranges from 0 (minimum dissimilarity) to 1.  It can be interpreted as the smallest fraction of units that need to be reclassified in order to make the distributions equal.  When \code{p2} is the reference distribution (true or expected distribution under a given hypothesis) than, following the Agresti's rule of thumb (Agresti 2002, pp. 329--330) , values of \eqn{\Delta_{12}<0.03}{D < 0.03} denotes that the estimated distribution \code{p1} follow the true or expected pattern quite closely.

\emph{Overlap} between two distributions:

\deqn{O_{12} = \sum_{j=1}^J min(p_{1,j},p_{2,j})  }{O = sum_j min(p_1,j , p_2,j) }

It is a measure of similarity which ranges from 0 to 1 (the distributions are equal).  It is worth noting that \eqn{O_{12}=1-\Delta_{12}}{O = 1 - D}.


\emph{Bhattacharyya coefficient}:

\deqn{B_{12} = \sum_{j=1}^J \sqrt{p_{1,j} \times p_{2,j}}  }{B = sum_j sqrt(p_1,j * p_2,j) }

It is a measure of similarity and ranges from 0 to 1 (the distributions are equal).

\emph{Hellinger's distance}:

\deqn{d_{H,12} = \sqrt{1-B_{1,2}} }{d_H = sqrt(1 - B) }

It is a dissimilarity measure which ranges from 0 (distributions are equal) to 1 (max dissimilarity).  It satisfies all the properties of a distance measure (\eqn{0 \leq d_{H,12} \leq 1}{0 <= d_H <= 1};  symmetry and triangle inequality).  
Hellinger's distance is related to the dissimilarity index, and it is possible to show that:

\deqn{d_{H,12}^2 \leq  \Delta_{12} \leq d_{H,12}\sqrt{2} }{ d_H^2 <= D <= d_H * sqrt(2) }


Alongside with those similarity/dissimilarity measures the Pearson's Chi-squared is computed. Two formulas are considered.  When \code{p2} is the reference distribution (true or expected under some hypothesis, \code{ref=TRUE}):

\deqn{ \chi^2_P = n_1 \sum_{j=1}^J \frac{\left( p_1,j - p_{2,j}\right)^2}{p_{2,j}}  }{ Chi_P = n_1 * sum_j (p_1,j - p_2,j)^2/(p_2,j) }

When \code{p2} is a distribution estimated on a second sample then:

\deqn{ \chi^2_P = \sum_{i=1}^2 \sum_{j=1}^J n_i \frac{\left( p_{i,j} - p_{+,j}\right)^2}{p_{+,j}}  }{ Chi_P = sum_i * sum_j n_i * (p_i,j - p_+,j)^2/(p_+,j) }

where \eqn{p_{+,j}}{p_+,j} is the expected frequency for category j, obtained as follows:

\deqn{ p_{+,j} = \frac{n_1 p_{1,j} + n_2 p_{2,j}}{n_1+n_2}  }{ p_+,j = (n_1*p_1,j + n_2*p_2,j)/(n_1+n_2) }

The Chi-Square value can be used to test the hypothesis that two distributions are equal (df=J-1).  Unfortunately such a test would not be useful when the distribution are estimated from samples selected from finite population using complex selection schemes (stratification, clustering, etc.).  In such a case different alternative corrected Chi-square tests are available (cf. Sarndal et al., 1992, Sec. 13.5). One possibility consist in dividing the Pearson's Chi-square test by the \emph{generalised design effect} of both the surveys.  Its estimation is not simple (sampling design variables need to be available). The generalised design effect is smaller than 1 in the presence of stratified random sampling designs.  It exceeds 1 the presence of a two stage cluster sampling design.  For the purposes of analysis it is reported the value of the generalised design effect \emph{g} that would determine the acceptance of the null hypothesis (equality of distributions) in the case of alpha=0.05 (df=J-1), i.e. values of \emph{g} such that 

\deqn{   \frac{\chi^2_P}{g} \leq \chi^2_{J-1,0.05} }{ Chi_P/d <= Chi_(J-1,0.05) }

}

\value{
A \code{list} object with two or three components depending on the argument \code{ref}. 

\item{meas}{A vector with the measures of similarity/dissimilarity between the distributions: dissimilarity index (\code{"tvd"}), overlap (\code{"overlap"}), Bhattacharyya coefficient \cr 
(\code{"Bhatt"}) and Hellinger's distance (\code{"Hell"}).}

\item{chi.sq}{ A vector with the following values: Pearson's Chi-square (\code{"Pearson"}), the degrees of freedom (\code{"df"}), the percentile of a Chi-squared distribution (\code{"q0.05"}) and the largest admissible value of the generalised design effect that would determine the acceptance of H0 (equality of distributions).}

\item{p.exp}{ When \code{ref=FALSE} it is reported the value of the reference distribution \eqn{p_{+,j}}{p_+,j} estimated used in deriving the Chi-square statistic and also the dissimilarity index.  On the contrary (\code{ref=FALSE})  it is set equal to the argument \code{p2}.}
}

\references{
Agresti A (2002) \emph{Categorical Data Analysis. Second Edition}. Wiley, new York.

Sarndal CE, Swensson B, Wretman JH (1992) \emph{Model Assisted Survey Sampling}. Springer--Verlag, New York.
}

\author{
 Marcello D'Orazio \email{madorazi@istat.it} 
}

%\seealso{ 
%
%}

\examples{
data(quine, package="MASS") #loads quine from MASS
str(quine)

# split quine in two subsets
set.seed(124)
lab.A <- sample(nrow(quine), 70, replace=TRUE)
quine.A <- quine[lab.A, c("Eth","Sex","Age")]
quine.B <- quine[-lab.A, c("Eth","Sex","Age")]

# compare est. distributions from 2 samples
# 1 variable
tt.A <- xtabs(~Age, data=quine.A)
tt.B <- xtabs(~Age, data=quine.B)
comp.prop(p1=tt.A, p2=tt.B, n1=nrow(quine.A), n2=nrow(quine.B), ref=FALSE)

# joint distr. of more variables
tt.A <- xtabs(~Eth+Sex+Age, data=quine.A)
tt.B <- xtabs(~Eth+Sex+Age, data=quine.B)
comp.prop(p1=tt.A, p2=tt.B, n1=nrow(quine.A), n2=nrow(quine.B), ref=FALSE)

# compare est. distr. with a one considered as reference
tt.A <- xtabs(~Eth+Sex+Age, data=quine.A)
tt.all <- xtabs(~Eth+Sex+Age, data=quine)
comp.prop(p1=tt.A, p2=tt.all, n1=nrow(quine.A), n2=NULL, ref=TRUE)


}

\keyword{multivariate}
