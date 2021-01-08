% --- Source file: RANDwNND_hotdeck.Rd ---
\name{RANDwNND.hotdeck}
\alias{RANDwNND.hotdeck}
\title{Random Distance hot deck.}

\description{
  This function implements a variant of the distance hot deck method.  For each recipient record a subset of of the closest donors is retained and then a donor is selected at random.
}

\usage{
RANDwNND.hotdeck(data.rec, data.don, match.vars=NULL, 
                 don.class=NULL, dist.fun="Manhattan", 
                 cut.don="rot", k=NULL, weight.don=NULL, 
                 keep.t=FALSE, ...)
}

\arguments{

\item{data.rec}{
  	A numeric matrix or data frame that plays the role of \emph{recipient}.  This data frame must contain the variables (columns), specified via \code{match.vars} and \code{don.class}, that should be used in the matching. 
  	
Missing values (\code{NA}) are allowed.
}

\item{data.don}{
   A matrix or data frame that plays the role of \emph{donor}.  This data frame must contain the variables (columns), specified via \code{match.vars} and \code{don.class}, that should be used in the matching. 
}

\item{match.vars}{
A character vector with the names of the variables (the columns in both the data frames) that have to be used to compute distances between records (rows) in \code{data.rec} and those in \code{data.don}.  When no matching variables are considered (\code{match.vars=NULL}) then all the units in the same donation class are considered as possible donors.  Hence one of them is selected at random or with probability proportional to its weight (see argument \code{weight.don}).  When \code{match.vars=NULL} and the donation classes are not created \cr
(\code{don.class=NULL}) then all the available records in the \code{data.don} are considered as potential donors.
}

\item{don.class}{
A character vector with the names of the variables (columns in both the data frames) that have to be used to identify donation classes.  In this case the computation of distances is limited to those units in \code{data.rec} and \code{data.doc} that belong to the same donation class.  The case of empty donation classes should be avoided.  It would be preferable that variables used to form donation classes are defined as \code{factor}.

When not specified (default), no donation classes are used.  This may result in a heavy computational effort.
}

\item{dist.fun}{
A string with the name of the distance function that has to be used.  The following distances can be used: \dQuote{Manhattan} (aka \dQuote{City block}; default), \dQuote{Euclidean}, \dQuote{Mahalanobis},\dQuote{exact} or \dQuote{exact matching}, \dQuote{Gower}, \dQuote{minimax}, \dQuote{difference}, or one of the distance functions available in the package \pkg{proxy}.  Note that the distances are computed using the function  \code{\link[proxy]{dist}} of the package \pkg{proxy} with the exception of the \dQuote{Gower} (see function \code{\link[StatMatch]{gower.dist}} for details), \dQuote{Mahalanobis} (function \code{\link[StatMatch]{mahalanobis.dist}}), \dQuote{minimax} (see \code{\link[StatMatch]{maximum.dist}}) \dQuote{difference} case.  Note that \code{dist.fun="difference"} computes just the difference between the values of the unique numeric matching variable considered; in practice, it should be used when the subset of the donation classes should be formed by comparing the values of the unique matching variable (for further details see the argument \code{cut.don}.

By setting \code{dist.fun="ANN"} or \code{dist.fun="RANN"} it is possible to search for the \code{k} nearest neighbours for each recipient record by using the the Approximate Nearest Neighbor (ANN) search as implemented in the  function \code{\link[RANN]{nn2}} provided by the package \pkg{RANN}.

When \code{dist.fun="Manhattan"}, \code{"Euclidean"}, \code{"Mahalanobis"} or \code{"minimax"} all the variables in \code{data.rec} and \code{data.don} must be numeric.  On the contrary, when \code{dist.fun="exact"} or \cr
\code{dist.fun="exact matching"}, all the variables in \code{data.rec} and \code{data.don} will be converted to character and, as far as the distance computation is concerned, they will be treated as categorical nominal variables, i.e. distance is 0 if a couple of units shows the same response category and 1 otherwise. 
}

\item{cut.don}{
A character string that, jointly with the argument \code{k}, identifies the rule to be used to form the subset of the closest donor records. 
\itemize{

			\item \code{cut.don="rot"}: (default) then the number of the closest donors to retain is given by \eqn{ \left[ \sqrt{n_{D}} \right]+1}{([n_D+1])^(1/2)}; being \eqn{ n_{D} }{n_D} the total number of available donors.  In this case \code{k} must not to be specified.
			
			\item \code{cut.don="span"}: the number of closest donors is determined as the proportion \code{k} of all the available donors, i.e. \eqn{ \left[ n_{D} \times k \right] }{[n_D * k]}.  Note that, in this case, \eqn{ 0< \code{k} \leq 1 }{0 < k <= 1}.
			
			\item \code{cut.don="exact"}: the \code{k}th closest donors out of the \eqn{n_{D}}{n_D} are retained.  In this case, \eqn{ 0< \code{k} \leq{ n_{D} } }{0 < k <= n_D}.
			
			\item \code{cut.don="min"}: the donors at the minimum distance from the recipient are retained.
			
			\item \code{cut.don="k.dist"}: only the donors whose distance from the recipient is less or equal to the value specified with the argument \code{k}.  Note that in this case it is not possible to use \code{dist.fun="ANN"}.
			
			\item \code{cut.don="lt"} or \code{cut.don="<"}: only the donors whose value of the matching variable is smaller than the value of the recipient are retained. Note that in this case it is has to be set \code{dist.fun="difference"}.
			
			\item \code{cut.don="le"} or \code{cut.don="<="}: only the donors whose value of the matching variable is smaller or equal to the value of the recipient are retained.  Note that in this case it is has to be set \code{dist.fun="difference"}.
			
			\item \code{cut.don="ge"} or \code{cut.don=">="}: only the donors whose value of the matching variable is greater or equal to the value of the recipient are retained.  Note that in this case it is has to be set \code{dist.fun="difference"}.
			
			\item \code{cut.don="gt"} or \code{cut.don=">"}: only the donors whose value of the matching variable is greater than the value of the recipient are retained. Note that in this case it is has to be set \code{dist.fun="difference"}.

}

}
\item{k}{
Depends on the \code{cut.don} argument. 
}

\item{weight.don}{ 
A character string providing the name of the variable with the weights associated to the  donor units in \code{data.don}.  When this variable is specified, then the selection of a donor among those in the subset of the closest donors is done with probability proportional to its weight (units with larger weight will have a higher chance of being selected).  When \code{weight.don=NULL} (default) all the units in the subset of the closest donors will have the same probability of being selected.
}

\item{keep.t}{
Logical, when donation classes are used by setting  \code{keep.t=TRUE} prints information on the donation classes being processed (by default \code{keep.t=FALSE}).
}

\item{...}{ 
Additional arguments that may be required by \code{\link[StatMatch]{gower.dist}}, by \cr
\code{\link[StatMatch]{maximum.dist}}, by \code{\link[proxy]{dist}} or by \code{\link[RANN]{nn2}}.
}

}
  
\details{
This function finds a donor record for each record in the recipient data set.  The donor is chosen at random in the subset of available donors.  This procedure is known as \emph{random hot deck} (cf. Andridge and Little, 2010).  In \code{RANDwNND.hotdeck}, the number of closest donors retained to form the subset is determined according to criterion specified with the argument \code{cut.don}.
The selection of the donor among those in the subset is carried out with equal probability (\code{weight.don=NULL}) or with probability proportional to a weight associated to the donors, specified via the \code{weight.don} argument.  This procedure is is known as \emph{weighted random hot deck} (cf. Andridge and Little, 2010).

The search for the subset of the closest donors can be speed up by using the Approximate Nearest Neighbor search as implemented in the  function \code{\link[RANN]{nn2}} provided by the package \pkg{RANN}.  Note that this search can be used in all the cases with the exception of \code{cut.don="k.dist"}.

Note that the same donor can be used more than once.

This function can also be used to impute missing values in a data set.  In this case \code{data.rec} is the part of the initial data set that contains missing values; on the contrary, \code{data.don} is the part of the data set without missing values. See \R code in the Examples for details.
} 


\value{

A \R list with the following components:

\item{mtc.ids}{
A matrix with the same number of rows of \code{data.rec} and two columns.  The first column contains the row names of the \code{data.rec} and the second column contains the row names of the corresponding donors selected from the \code{data.don}.  When the input matrices do not contain row names, then a numeric matrix with the indexes of the rows is provided.
}

\item{sum.dist}{
A matrix with summary statistics concerning the subset of the closest donors.  The first three columns report the minimum, the maximum and the standard deviation of the distances among the recipient record and the donors in the subset of the closest donors, respectively.  The 4th column reports the cutting distance, i.e. the value of the distance such that donors at a higher distance are discarded.  The 5th column reports the distance between the recipient and the donor chosen at random in the subset of the donors.
}

\item{noad}{
For each recipient unit, reports the number of donor records in the subset of closest donors. 
}

\item{call}{
How the function has been called.
}

}

\references{

Andridge, R.R., and Little, R.J.A. (2010) \dQuote{A Review of Hot Deck Imputation for Survey Non-response}. \emph{International Statistical Review}, \bold{78}, 40--64. 

D'Orazio, M., Di Zio, M. and Scanu, M. (2006). \emph{Statistical Matching: Theory and Practice.} Wiley, Chichester.

Rodgers, W.L. (1984). \dQuote{An evaluation of statistical matching}. \emph{Journal of Business and Economic Statistics}, \bold{2}, 91--102.

Singh, A.C., Mantel, H., Kinack, M. and Rowe, G. (1993). \dQuote{Statistical matching: use of auxiliary information as an alternative to the conditional independence assumption}. \emph{Survey Methodology}, \bold{19}, 59--79.
}

\author{
 Marcello D'Orazio \email{mdo.statmatch@gmail.com} 
}

\seealso{ 
 
\code{\link[StatMatch]{NND.hotdeck}}
}

\examples{

data(samp.A, samp.B, package="StatMatch") #loads data sets
?samp.A
?samp.B


# samp.A plays the role of recipient
# samp.B plays the role of donor
# find a donor in the in the same region ("area5") and with the same
# gender ("sex"), then only the closest k=20 donors in terms of 
# "age" are cnsidered and one of them is picked up at random

out.RND.1 <- RANDwNND.hotdeck(data.rec=samp.A, data.don=samp.B,
                              don.class=c("area5", "sex"), dist.fun="ANN",
                              match.vars="age", cut.don="exact", k=20)

# create the synthetic (or fused) data.frame:
# fill in "labour5" in A
fused.1 <- create.fused(data.rec=samp.A, data.don=samp.B,
                        mtc.ids=out.RND.1$mtc.ids, z.vars="labour5")
head(fused.1)

# weights ("ww") are used in selecting the donor in the final step

out.RND.2 <- RANDwNND.hotdeck(data.rec=samp.A, data.don=samp.B,
                              don.class=c("area5", "sex"), dist.fun="ANN",
                              match.vars="age", cut.don="exact", 
                              k=20, weight.don="ww")
fused.2 <- create.fused(data.rec=samp.A, data.don=samp.B,
                        mtc.ids=out.RND.2$mtc.ids, z.vars="labour5")
head(fused.2)

# find a donor in the in the same region ("area5") and with the same
# gender ("sex"), then only the donors with "age" <= to the age of the
# recipient are considered,
# then one of them is picked up at random

out.RND.3 <- RANDwNND.hotdeck(data.rec=samp.A, data.don=samp.B,
                              don.class=c("area5", "sex"), dist.fun="diff",
                              match.vars="age", cut.don="<=")

# create the synthetic (or fused) data.frame:
# fill in "labour5" in A
fused.3 <- create.fused(data.rec=samp.A, data.don=samp.B,
                        mtc.ids=out.RND.3$mtc.ids, z.vars="labour5")
head(fused.3)

# Example of Imputation of missing values
# introducing missing vales in iris
ir.mat <- iris
miss <- rbinom(nrow(iris), 1, 0.3)
ir.mat[miss==1,"Sepal.Length"] <- NA
iris.rec <- ir.mat[miss==1,-1]
iris.don <- ir.mat[miss==0,]

#search for NND donors
imp.RND <- RANDwNND.hotdeck(data.rec=iris.rec, data.don=iris.don,
                            match.vars=c("Sepal.Width","Petal.Length", "Petal.Width"),
                            don.class="Species")

# imputing missing values
iris.rec.imp <- create.fused(data.rec=iris.rec, data.don=iris.don,
                             mtc.ids=imp.RND$mtc.ids, z.vars="Sepal.Length")

# rebuild the imputed data.frame
final <- rbind(iris.rec.imp, iris.don)
head(final)

}

\keyword{nonparametric}