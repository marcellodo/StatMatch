% --- Source file: RANDwNND_hotdeck.Rd ---
\name{RANDwNND.hotdeck}
\alias{RANDwNND.hotdeck}
\title{Random Distance hot deck method for Statistical Matching.}

\description{
  This function implements a variant of the distance hot deck method. For each recipient record a subset of units consisting of the closest donors is retained and then a donor is selected at random
}

\usage{
RANDwNND.hotdeck(data.rec, data.don, match.vars, don.class=NULL,  
                   dist.fun="Euclidean", cut.don="rot", k=NULL) 
}

\arguments{

\item{data.rec}{
  	A numeric matrix or data frame that plays the role of \emph{recipient} in the statistical matching application. This data frame must contain the variables (columns) that should be used, directly or indirectly, in the computation of distances between its observations (rows) and those in \code{data.y}. 
  	
Missing values (\code{NA}) are allowed.
}

\item{data.don}{
   A matrix or data frame that plays the role of \emph{donor} in the statistical matching application. The variables (columns) involved, directly or indirectly, in the computation of distances be the same and of the same type as those in \code{data.rec}.   
}

\item{match.vars}{
A character vector with the names of the variables (the columns in both the data frames) that have to be used to compute distances among records (rows) in \code{data.rec} and those in \code{data.don}.
}

\item{don.class}{
A character vector with the names of the variables (columns in both the data frames) that have to be used to identify donation classes. In this case the computation of distances is limited to those units in \code{data.rec} and \code{data.doc} that belong to the same donation classes. The case of empty donation classes should be avoided. To avoid confusion, it is preferable the variables used to form donation classes are defined as \code{factor}.

When not specified (default), no donation classes are used. This may result in a heavy computational effort.
}

\item{dist.fun}{
A string with the name of the distance function that has to be used. The following distances are allowed: \dQuote{Euclidean} (default), \dQuote{Manhattan} (aka \dQuote{City block}), \dQuote{exact} or \dQuote{exact matching}, \dQuote{Gower} or one of the distance functions available in the package \pkg{proxy}. Note that the distance is computed using the function  \code{\link[proxy]{dist}} of the package \pkg{proxy} with the exception of the \dQuote{Gower} case. 

When \code{dist.fun= "Euclidean" } or \dQuote{Manhattan} all the non numeric variables in \code{data.rec} and \code{data.don} will be converted to numeric. On the contrary, when \code{dist.fun="exact"} or \dQuote{exact matching}, all the variables in \code{data.rec} and \code{data.don} will be converted to character and, as far as distance computation is concerned, they will be treated as categorical nominal variables, i.e.\ the distance is 0 if two units show the same response category and 1 otherwise. 

When \code{dist.fun="Gower"} the Gower dissimilarity is considered. See function \code{\link[StatMatch]{gower.dist}} for details.
}

\item{cut.don}{
A character string that, jointly with the argument \code{k}, identifies the rule to be used to form the subset of the closest donor records. 
\itemize{

			\item \code{cut.don="rot"}: (default) then the number of the closest donors to retain is given by \eqn{ \left[ \sqrt{n_{D}} \right]+1}{([n_D+1])^(1/2)}; being \eqn{ n_{D} }{n_D} the total number of available donors. In this case \code{k} must not to be specified.
			
			\item \code{cut.don="span"}: the number of closest donors is determined as the proportion \code{k} of all the available donors, i.e. \eqn{ \left[ n_{D} \times k \right] }{[n_D * k]}. Note that, in this case, \eqn{ 0< \code{k} \leq 1 }{0 < k <= 1}.
			
			\item \code{cut.don="exact"}: the \code{k}th closest donors out of the \eqn{n_{D}}{n_D} are retained. In this case, \eqn{ 0< \code{k} \leq{ n_{D} } }{0 < k <= n_D}.
			
			\item \code{cut.don="min"}: the donors at the minimum distance from the recipient are retained.
			
			\item \code{cut.don="k.dist"}: only the donors whose distance from the recipient is less or equal to the value specified with the argument \code{k}.
			
}

}
\item{k}{
Depends on the \code{cut.don} argument. 
}
}
  
\details{
This function finds a donor record for each record in the recipient data set. This donor is chosen at random in the subset of available donors. The number of closest donors retained to form the subset is determined according to criterion specified with the argument \code{cut.don}.

Note that the same donor can be used more than once.

Note that this function van also be used to impute missing values in a data set. In this case \code{data.rec} is the part of the initial data set that contains missing values; on the contrary, \code{data.don} is the part of the data set without missing values. See \R code in the Examples for details.
} 


\value{

A \R list with the following components:

\item{mtc.ids}{
A matrix with the same number of rows of \code{data.rec} and two columns. The first column contains the row names of the \code{data.rec} and the second column contains the row names of the corresponding donors selected from the \code{data.don}. When the input matrices do not contain row names, then a numeric matrix with the indexes of the rows is provided.
}

\item{sum.dist}{
A matrix with summary statistics concerning the subset of the closest donors from which a donor for a given recipient is chosen. The first three columns report the minimum, the maximum and the standard deviation of the distances among the recipient record and the donors in the subset of the closest donors, respectively. The 4th column reports the cutting distance, i.e. the value of the distance such that donors at a higher distance are discarded. The 5th column reports the distance between the recipient and the donor chosen at random in the subset of the donors.
}

\item{noad}{
For each recipient unit, reports the number of donor records in the subset of closest donors. 
}

\item{call}{
How the function has been called.
}

}

\references{

D'Orazio, M., Di Zio, M. and Scanu, M. (2006). \emph{Statistical Matching: Theory and Practice.} Wiley, Chichester.

Rodgers, W.L. (1984). \dQuote{An evaluation of statistical matching}. \emph{Journal of Business and Economic Statistics}, \bold{2}, 91--102.

Singh, A.C., Mantel, H., Kinack, M. and Rowe, G. (1993). \dQuote{Statistical matching: use of auxiliary information as an alternative to the conditional independence assumption}. \emph{Survey Methodology}, \bold{19}, 59--79.
}

\author{
 Marcello D'Orazio \email{madorazi@istat.it} 
}

\seealso{ 
 
\code{\link[StatMatch]{NND.hotdeck}}
}

\examples{

# reproduce the classical matching framework
lab <- c(1:10, 51:60, 101:110)
iris.rec <- iris[lab, c(1:3,5)]  # recipient data.frame
iris.don <- iris[-lab, c(1:2,4:5)] # recipient data.frame

# Now iris.rec and iris.don have the variables
# "Sepal.Length", "Sepal.Width"  and "Species"
# in common.
#  "Petal.Length" is available only in iris.rec
#  "Petal.Width"  is available only in iris.don

# find a donor in the subset of closest donors using cut.don="rot";
# the distance is computed using "Sepal.Length" and "Sepal.Width"

out.NND.1 <- RANDwNND.hotdeck(data.rec=iris.rec, data.don=iris.don,
              match.vars=c("Sepal.Length", "Sepal.Width") )

# create the synthetic (or fused) data.frame:
# fill in "Petal.Width" in iris.rec
fused.1 <- create.fused(data.rec=iris.rec, data.don=iris.don, 
             mtc.ids=out.NND.1$mtc.ids, z.vars="Petal.Width") 


# find a donor in the subset of closest donors using cut.don="rot";
# the distance is computed using "Sepal.Length" and "Sepal.Width"
# "Species" is used to form donation classes

out.NND.2 <- RANDwNND.hotdeck(data.rec=iris.rec, data.don=iris.don,
              match.vars=c("Sepal.Length", "Sepal.Width") , don.class="Species")


# as before, but with a different criteria to reduce the no. of donors:
# the first half (k=0.5) of the closest available donors is retained,
# then a donor is chosen at random among them

out.NND.3 <- RANDwNND.hotdeck(data.rec=iris.rec, data.don=iris.don,
              don.class="Species", match.vars=c("Sepal.Length", "Sepal.Width"),
              cut.don="span", k=0.5 )


# as before, but the subset of closest donors is formed by considering
# only the first k=5 closest donors

out.NND.4 <- RANDwNND.hotdeck(data.rec=iris.rec, data.don=iris.don,
              don.class="Species", match.vars=c("Sepal.Length", "Sepal.Width"),
              cut.don="exact", k=5 )


# as before, but the subset of closest donors is formed by considering
# the donors whose distance (Gower) is less or equal to k=0.33

out.NND.5 <- RANDwNND.hotdeck(data.rec=iris.rec, data.don=iris.don,
              don.class="Species", match.vars=c("Sepal.Length", "Sepal.Width"),
              dist.fun="Gower", cut.don="k.dist", k=0.33 )


# Example of Imputation of missing values
# introducing missing vales in iris
ir.mat <- iris
miss <- rbinom(nrow(iris), 1, 0.3)
ir.mat[miss==1,"Sepal.Length"] <- NA
iris.rec <- ir.mat[miss==1,-1]
iris.don <- ir.mat[miss==0,]

#search for NND donors
imp.NND <- RANDwNND.hotdeck(data.rec=iris.rec, data.don=iris.don,
               match.vars=c("Sepal.Width","Petal.Length", "Petal.Width"),
               don.class="Species")

# imputing missing values
iris.rec.imp <- create.fused(data.rec=iris.rec, data.don=iris.don, 
             mtc.ids=imp.NND$mtc.ids, z.vars="Sepal.Length") 

# rebuild the imputed data.frame
final <- rbind(iris.rec.imp, iris.don)

}

\keyword{nonparametric}