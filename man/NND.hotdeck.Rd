% --- Source file: NND_hotdeck.Rd ---
\name{NND.hotdeck}
\alias{NND.hotdeck}
\title{Distance Hot Deck method.}

\description{
  This function implements the distance hot deck method to match the records of two data sources that share some variables. 
}

\usage{
NND.hotdeck(data.rec, data.don, match.vars, 
             don.class=NULL, dist.fun="Manhattan",
             constrained=FALSE, constr.alg="Hungarian", ...) 
}

\arguments{

\item{data.rec}{
  	A matrix or data frame that plays the role of \emph{recipient}.  This data frame must contain the variables (columns) that should be used, directly or indirectly, in the matching application (specified via \code{match.vars} and eventually \code{don.class}). 

Missing values (\code{NA}) are allowed.
}

\item{data.don}{
   A matrix or data frame that plays the role of \emph{donor}.  The variables (columns) involved, directly or indirectly, in the computation of distance must be the same and of the same type as those in \code{data.rec} (specified via \code{match.vars} and eventually \code{don.class}).   
}

\item{match.vars}{
A character vector with the names of the matching variables (the columns in both the data frames) that have to be used to compute distances among records (rows) in \code{data.rec} and those in \code{data.don}.
}

\item{don.class}{
A character vector with the names of the variables (columns in both the data frames) that have to be used to identify the donation classes.  In this case the computation of distances is limited to those units of \code{data.rec} and \code{data.doc} that belong to the same donation class.  The case of empty donation classes should be avoided.  It would be preferable that variables used to form donation classes are defined as \code{factor}.

When not specified (default) no donation classes are used.  This should require more memory to store the distance matrix and a higher computational effort.
}

\item{dist.fun}{
A string with the name of the distance function that has to be used.  The following distances are allowed: \dQuote{Manhattan} (aka \dQuote{City block}; default), \dQuote{Euclidean}, \dQuote{Mahalanobis},\dQuote{exact} or \dQuote{exact matching}, \dQuote{Gower}, \dQuote{minimax} or one of the distance functions available in the package \pkg{proxy}.  Note that the distance is computed using the function  \code{\link[proxy]{dist}} of the package \pkg{proxy} with the exception of the \dQuote{Gower} (see function \code{\link[StatMatch]{gower.dist}} for details), \dQuote{Mahalanobis} (function \code{\link[StatMatch]{mahalanobis.dist}}) and \dQuote{minimax} (see \code{\link[StatMatch]{maximum.dist}}) cases.

When \code{dist.fun="Manhattan"}, \code{"Euclidean"}, \code{"Mahalanobis"} or \code{"minimax"} all the non numeric variables in \code{data.rec} and \code{data.don} will be converted to numeric.  On the contrary, when \code{dist.fun="exact"} or \code{dist.fun="exact matching"}, all the variables in \code{data.rec} and \code{data.don} will be converted to character and, as far as the distance computation is concerned, they will be treated as categorical nominal variables, i.e. distance is 0 if a couple of units present the same response category and 1 otherwise. 

}

\item{constrained}{
Logical.  When \code{constrained=FALSE} (default) each record in \code{data.don} can be used as a donor more than once.  On the contrary, when \cr \code{constrained=TRUE} each record in \code{data.don} can be used as a donor only once.  In this case, the set of donors is selected by solving a transportation problem, in order to minimize the overall matching distance.  See description of the argument \code{constr.alg} for details.
}

\item{constr.alg}{
A string that has to be specified when \code{constrained=TRUE}.  Two choices are available: \dQuote{lpSolve} and \dQuote{hungarian}.  In the first case, \code{constr.alg="lpSolve"}, the transportation problem is solved by means of the function \code{\link[lpSolve]{lp.transport}} available in the package \pkg{lpSolve}.  When \code{constr.alg="hungarian"} (default) the transportation problem is solved using the Hungarian method, implemented in function \code{\link[clue]{solve_LSAP}} available in the package \pkg{clue}.  Note that Hungarian algorithm is faster and more efficient if compared to \code{constr.alg="lpSolve"} .
}

\item{...}{ 
Additional arguments that may be required by \code{\link[StatMatch]{gower.dist}}, or by \cr 
\code{\link[StatMatch]{maximum.dist}}, or by \code{\link[proxy]{dist}}.
}
}
  
\details{
This function finds a donor record in \code{data.don} for each record in \code{data.rec}.  In the unconstrained case, it searches for the closest donor record in \code{data.don} for each record in the recipient data set, according to the chosen distance function.  When for a given recipient record there are more donors available at the minimum distance, one of them is picked at random.

In the constrained case the set of donors is chosen in order to minimize the overall matching distance.  In this case the number of units (rows) in the donor data set has to be larger or equal to the number of units of the recipient data set.  When the donation classes are used, this condition must be satisfied in each donation class.  For further details on nearest neighbor distance hot deck refer to Chapter 2 in D'Orazio \emph{et al.} (2006).

This function can also be used to impute missing values in a data set using the nearest neighbor distance hot deck.  In this case \code{data.rec} is the part of the initial data set that contains missing values; on the contrary, \code{data.don} is the part of the data set without missing values.  See \R code in the Examples for details.
} 

\value{

A \R list with the following components:

\item{mtc.ids}{
A matrix with the same number of rows of \code{data.rec} and two columns.  The first column contains the row names of the \code{data.rec} and the second column contains the row names of the corresponding donors selected from the \code{data.don}.  When the input matrices do not contain row names, a numeric matrix with the indexes of the rows is provided.
}

\item{dist.rd}{
A vector with the distances among each recipient unit and the corresponding donor. 
}

\item{noad}{
When \code{constrained=FALSE}, it reports the number of available donors at the minimum distance for each recipient unit. 
}

\item{call}{
How the function has been called.
}

}

\references{

D'Orazio, M., Di Zio, M. and Scanu, M. (2006). \emph{Statistical Matching: Theory and Practice.} Wiley, Chichester.

Hornik K. (2012).  clue: Cluster ensembles.  R package version 0.3-45.  \url{http://CRAN.R-project.org/package=clue}.

Rodgers, W.L. (1984). \dQuote{An evaluation of statistical matching}. \emph{Journal of Business and Economic Statistics}, \bold{2}, 91--102.

Singh, A.C., Mantel, H., Kinack, M. and Rowe, G. (1993). \dQuote{Statistical matching: use of auxiliary information as an alternative to the conditional independence assumption}. \emph{Survey Methodology}, \bold{19}, 59--79.
}

\author{
 Marcello D'Orazio \email{madorazi@istat.it} 
}

\seealso{ 
\code{\link[StatMatch]{RANDwNND.hotdeck}}
}

\examples{

# reproduce the classical matching framework
lab <- c(1:15, 51:65, 101:115)
iris.rec <- iris[lab, c(1:3,5)]  # recipient data.frame 
iris.don <- iris[-lab, c(1:2,4:5)] #donor data.frame

# Now iris.rec and iris.don have the variables
# "Sepal.Length", "Sepal.Width"  and "Species"
# in common.
#  "Petal.Length" is available only in iris.rec
#  "Petal.Width"  is available only in iris.don

# Find the closest donors donors computing distance
# on "Sepal.Length" and "Sepal.Width"
# unconstrained case, Euclidean distance

out.NND.1 <- NND.hotdeck(data.rec=iris.rec, data.don=iris.don,
                match.vars=c("Sepal.Length", "Sepal.Width") )

# create the synthetic data.set:
# fill in "Petal.Width" in iris.rec

fused.1 <- create.fused(data.rec=iris.rec, data.don=iris.don, 
             mtc.ids=out.NND.1$mtc.ids, z.vars="Petal.Width") 


# Find the closest donors computing distance
# on "Sepal.Length", "Sepal.Width" and Species;
# unconstrained case, Gower's distance

out.NND.2 <- NND.hotdeck(data.rec=iris.rec, data.don=iris.don,
              match.vars=c("Sepal.Length", "Sepal.Width", "Species"), 
              dist.fun="Gower")


# find the closest donors using "Species" to form donation classes
# and "Sepal.Length" and "Sepal.Width" to compute distance;
# unconstrained case.

out.NND.3 <- NND.hotdeck(data.rec=iris.rec, data.don=iris.don,
              match.vars=c("Sepal.Length", "Sepal.Width"),
              don.class="Species")


# find the donors using "Species" to form donation classes
# and "Sepal.Length" and "Sepal.Width" to compute distance;
# constrained case, "Hungarian" algorithm

library(clue)
out.NND.4 <- NND.hotdeck(data.rec=iris.rec, data.don=iris.don,
                         match.vars=c("Sepal.Length", "Sepal.Width"),
                         don.class="Species", constrained=TRUE, 
                         constr.alg="Hungarian")


# find the donors using "Species" to form donation classes
# and "Sepal.Length" and "Sepal.Width" to compute distance;
# constrained case, transportation problem solved by functions 
# in package "lpSolve"

library(lpSolve)
out.NND.5 <- NND.hotdeck(data.rec=iris.rec, data.don=iris.don,
                         match.vars=c("Sepal.Length", "Sepal.Width"),
                         don.class="Species", constrained=TRUE, 
                         constr.alg="lpSolve")


# Example of Imputation of missing values.
# Introducing missing values in iris
ir.mat <- iris
miss <- rbinom(nrow(iris), 1, 0.3)
ir.mat[miss==1,"Sepal.Length"] <- NA
iris.rec <- ir.mat[miss==1,-1]
iris.don <- ir.mat[miss==0,]

#search for NND donors
imp.NND <- NND.hotdeck(data.rec=iris.rec, data.don=iris.don,
               match.vars=c("Sepal.Width","Petal.Length", "Petal.Width"),
               don.class="Species")

# imputing missing values
iris.rec.imp <- create.fused(data.rec=iris.rec, data.don=iris.don, 
             mtc.ids=imp.NND$mtc.ids, z.vars="Sepal.Length") 

# rebuild the imputed data.frame
final <- rbind(iris.rec.imp, iris.don)

}
\keyword{nonparametric}