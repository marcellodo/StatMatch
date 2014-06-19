% --- Source file:  ---
\name{samp.B}
\docType{data}
\alias{samp.B}

\title{Artificial data set resembling EU--SILC survey}

\description{
This data set provides a limited number of variables observed at persons levels among those usually collected in the European Union Statistics on Income and Living Conditions Survey (EU--SILC). The data are artificially generated, just to show the application of the statistical matching techniques implemented in \pkg{StatMatch}.
}

\usage{
data(samp.B)
}

\format{
A data frame with 6686 observations and the following variables:
\describe{
	\item{HH.P.id}{unique unit identifier of the type \code{aa.bb} where \code{aa} identifies the Household while \code{bb} identifies the household member}
	\item{area5}{large geographic area, factor with 5 levels: `NE'=North--East, `NO'=North--West,
	`C'=center, `S'=South, `I'=islands}
	\item{urb}{Degree of urbanization, factor with 3 levels: `1'=densely populated area, `2'=intermediate area,
	`3'=thinly populated area}
	\item{hsize}{integer, size of the household in which the person lives}
	\item{hsize5}{factor with 5 levels derived from \code{hsize}, where the 5th level `>=5' denotes 5 and more people in the household}
	\item{age}{integer, the person's age}
	\item{c.age}{factor, age categorized in 5 classes}
	\item{sex}{factor, the person's gender: `1'=male, `2'=female}
	\item{marital}{factor, the person's marital status: `1'=never married, `2'=married, `3'=other (separated, widowed, divorced)}
	\item{edu7}{factor, the person's highest education level attained, follows the ISCED-97 categories:
	 `0'=pre--primary education, `1'=primary education, `2'=lower secondary education, `3'= (upper) secondary education, `4'= post--secondary non tertiary education, `5'=first stage of tertiary education (not leading directly to an advanced research qualification), `6'=second stage of tertiary education (leading to an advanced research qualification)}
	 \item{labour5}{the person's self--defined economic status, factor with 5 levels: `1'=employee working full--time or part--time, `2'=self--employed working full--time or part--time, `3'=unemployed, `4'=In retirement or in early retirement or has given up business, `5'=other status (student, permanent disabled, in compulsory military service, fulfilling domestic tasks, etc.)}
	 \item{ww}{numeric, the unit's weight}
}
}

\details{
Please note that this data set is just for illustrative purposes. The unit's weight do not reflect the Italian population size. The variables included are derived starting from the those usually observed in the EU--SILC survey. 
}

\source{
This data set is artificially created starting from the EU--SILC survey structure. 
}

\references{
\url{http://epp.eurostat.ec.europa.eu/portal/page/portal/microdata/eu_silc}
}

\keyword{datasets}