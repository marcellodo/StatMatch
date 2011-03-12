% --- Source file: rankNND_hotdeck.Rd ---
\name{rankNND.hotdeck}
\alias{rankNND.hotdeck}
\title{ Rank distance hot deck method for Statistical Matching.}

\description{
  This function implements rank  hot deck distance method.  For each recipient record the closest donors is chosen by considering the distance among their empirical cumulative distribution function.
}

\usage{
rankNND.hotdeck(data.rec, data.don, var.rec, var.don=var.rec, 
                 don.class=NULL,  weight.rec=NULL, weight.don=NULL) 
}

\arguments{

\item{data.rec}{
  A numeric matrix or data frame that plays the role of \emph{recipient} in the statistical matching application.  This data frame must contain the variable \code{var.rec} to be used in computing the empirical cumulative distribution function and eventually the variables that should be used to identify the donation classes (see argument \code{don.class}) and the case weights (see argument \code{weight.rec}). 

Missing values (\code{NA}) are not allowed.
}

\item{data.don}{
   A matrix or data frame that plays the role of \emph{donor} in the statistical matching application.  This data frame must contain the variable \code{var.don} to be used in computing the empirical cumulative distribution function and eventually the variables that should be used to identify the donation classes (see argument \code{don.class}) and the case weights (see argument \code{weight.don}).  
}

\item{var.rec}{
A character vector with the name of the variable that should be ranked. 
}

\item{var.don}{
A character vector with the name of the variable that should be ranked.  If not specified, by default \code{var.don=var.rec}.
}

\item{don.class}{
A character vector with the names of the variables (columns in both the data frames) that have to be used to identify donation classes.  In this case the computation of ranks it is limited to those units in \code{data.rec} and \code{data.doc} that belong to the same donation classes.  The case of empty donation classes should be avoided.  It would be preferable the variables used to form donation classes are defined as \code{factor}.

When not specified (default), no donation classes are used.
}

\item{weight.rec}{
Eventual name of the variable in \code{data.rec} that provides the weights that should be used in computing the the empirical cumulative distribution function for \code{var.rec} (see Details).
}

\item{weight.don}{
Eventual name of the variable in \code{data.don} that provides the weights that should be used in computing the the empirical cumulative distribution function for \code{var.don} (see Details).
}

}
  
\details{
This function finds a donor record for each record in the recipient data set.  The chosen donor is the one at the closest distance in terms of empirical cumulative distribution (Singh et al., 1990).  In practice the distance is computed by considering the estimated empirical cumulative distribution for the reference variable (\code{var.rec} and \code{var.don}) in \code{data.rec} and \code{data.don}.  The empirical cumulative distribution function is estimated by:

\deqn{ \hat{F}(y_k) = \frac{1}{n} \sum_{i=1}^{n} I(y_i\leq y_k) }{F(y_k) = (1/n) sum_i I(y_i <= y_k)}

being \eqn{I()=1}{I()=1} if \eqn{y_i\leq y_k}{ y_i <= y_k} and 0 otherwise.

In the presence of weights the empirical cumulative distribution function is estimated by:

\deqn{ \hat{F}(y_k) =  \frac{\sum_{i=1}^{n} w_i I(y_i\leq y_k)}{\sum_{i=1}^{n} w_i} }{F(y_k) = (sum_i w_i*I(y_i <= y_k)) /(sum_i w_i)}

When there are more donors at the same distance then one is chosen at random.

Note that when the donation classes are introduced then empirical cumulative distribution function is estimated independently in each donation classes and the search of a recipient is restricted to donors in the same donation class.

A donor can be chosen more than once.
} 


\value{

A \R list with the following components:

\item{mtc.ids}{
A matrix with the same number of rows of \code{data.rec} and two columns.  The first column contains the row names of the \code{data.rec} and the second column contains the row names of the corresponding donors selected from the \code{data.don}.  When the input matrices do not contain row names, then a numeric matrix with the indexes of the rows is provided.
}
\item{dist.rd}{
A vector with the distances among each recipient unit and the corresponding donor. 
}

\item{noad}{
The number of available donors at the minimum distance for each recipient unit. 
}

\item{call}{
How the function has been called.
}

}

\references{

D'Orazio, M., Di Zio, M. and Scanu, M. (2006). \emph{Statistical Matching: Theory and Practice.} Wiley, Chichester.

Singh, A.C., Mantel, H., Kinack, M. and Rowe, G. (1993). \dQuote{Statistical matching: use of auxiliary information as an alternative to the conditional independence assumption}. \emph{Survey Methodology}, \bold{19}, 59--79.
}

\author{
 Marcello D'Orazio \email{madorazi@istat.it} 
}

\seealso{ 
 
\code{\link[StatMatch]{NND.hotdeck}}
}

\examples{
require(SDaA)
data(agpop, agsrs, agstrat, package="SDaA") #loads ag datasets from SDaA
str(agpop)
str(agsrs)
str(agstrat)

agsrs$w.srs <- nrow(agpop)/nrow(agsrs) # add weights

# adds region to agsrs
state.region <- data.frame(xtabs(weight~state+region, data=agstrat))
state.region <- subset(state.region, Freq>0)
agsrs <- merge(agsrs, state.region[,1:2], by="state", all.x=TRUE)

# simulate statistical matching framework
A <- agsrs[, c("region", "acres82", "acres87", "w.srs")]
B <- agstrat[, c("region", "acres82", "acres92", "weight")]

# simplest call to rankNND.hotdeck()
out.1 <- rankNND.hotdeck(data.rec=A, data.don=B, var.rec="acres82")
fused.1 <- create.fused(data.rec=A, data.don=B,
             mtc.ids=out.1$mtc.ids, z.vars="acres92")
head(fused.1)

#  call to rankNND.hotdeck() with usage of weights
out.2 <- rankNND.hotdeck(data.rec=A, data.don=B, var.rec="acres82",
            weight.rec="w.srs", weight.don="weight")
fused.2 <- create.fused(data.rec=A, data.don=B,
             mtc.ids=out.2$mtc.ids, z.vars="acres92")
head(fused.2)

#  call to rankNND.hotdeck() with usage of weights and don classes
out.3 <- rankNND.hotdeck(data.rec=A, data.don=B, var.rec="acres82",
            don.class="region", weight.rec="w.srs", weight.don="weight")
fused.3 <- create.fused(data.rec=A, data.don=B,
             mtc.ids=out.3$mtc.ids, z.vars="acres92")
head(fused.3)

}

\keyword{nonparametric}