% --- Source file: create_fused.Rd ---
\name{create.imputed}
\alias{create.imputed}
\title{Fills-in missing values in the recipient dataset with values observed on the donors units}

\description{
Imputes the missing values (NAs) in the recipient dataset with values observed on the donors units after search of donors with NND or random hotdeck.
}

\usage{
create.imputed(data.rec, data.don, mtc.ids)  
}

\arguments{

\item{data.rec}{
A matrix or data frame that has missing values. 
}

\item{data.don}{
A matrix or data frame that is used for donation (imputation).
}

\item{mtc.ids}{
A matrix with two columns.  Each row must contain the name or the index of the recipient record (row) in \code{data.don} and the name or the index of the corresponding donor record (row) in \code{data.don}.  Note that this type of matrix is returned by the functions \code{\link[StatMatch]{NND.hotdeck}}, \code{\link[StatMatch]{RANDwNND.hotdeck}}, \code{\link[StatMatch]{rankNND.hotdeck}}, and \code{\link[StatMatch]{mixed.mtc}}.
}

}




\details{
This function allows to fill in the missing values (NAs) in the recipient with values observed in the donor data set after the search of donors via NND or random hotdeck with available functions in the package, i.e. \code{\link[StatMatch]{NND.hotdeck}}, \code{\link[StatMatch]{RANDwNND.hotdeck}}, \code{\link[StatMatch]{rankNND.hotdeck}}, and \code{\link[StatMatch]{mixed.mtc}}.
When the same recorc in the recipient dataset presents 2 or more NAs they all will be replaced with values observed on the chosen donor for that unit; this corresponds to joint hotdeck imputation. 
} 


\value{

The data frame \code{data.rec} missing values (NAs) filled in.

}
\references{
D'Orazio, M., Di Zio, M. and Scanu, M. (2006). \emph{Statistical Matching: Theory and Practice.} Wiley, Chichester.
}

\author{
 Marcello D'Orazio \email{mdo.statmatch@gmail.com} 
}

\seealso{ 
\code{\link[StatMatch]{NND.hotdeck}}
\code{\link[StatMatch]{RANDwNND.hotdeck}}
\code{\link[StatMatch]{rankNND.hotdeck}}

}

\examples{

# introduce missing values
# in Petal.Length of iris dataset

set.seed(13579)
pos <- sample(x = 1:nrow(iris), size = 15, 
              replace = FALSE)


iris.rec <- iris[pos, ]  # recipient data.frame with missing values
iris.rec[, "Petal.Length"] <- NA

iris.don <- iris[-pos, ] # donor data.frame ALL observed

# find the closest donors using NND hot deck;
# distances are computed on "Petal.Width"
# donors only of the same Specie

out.NND <- NND.hotdeck(data.rec=iris.rec, data.don=iris.don,
                       match.vars=c("Petal.Width"), 
                       don.class="Species")

# impute missing
iris.rec.imp <- create.imputed(data.rec=iris.rec, data.don=iris.don, 
                        mtc.ids=out.NND$mtc.ids)
summary(iris.rec.imp$Petal.Length)

}
\keyword{manip}