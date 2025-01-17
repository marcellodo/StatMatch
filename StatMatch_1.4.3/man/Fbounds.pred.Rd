% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Fbounds.pred.R
\name{Fbounds.pred}
\alias{Fbounds.pred}
\title{
Estimates Frechet bounds for cells in the contingency table crossing two categorical variables observed in distinct samples referred to the same target population.}
\description{
This function assesses the uncertainty in estimating the contingency table crossing \code{y.rec} (Y) and \code{z.don} (Z) when the two variables are observed in two different samples sharing a number of common predictors.}
\usage{
Fbounds.pred(data.rec, data.don,
             match.vars, y.rec, z.don, pred = "multinom",
             w.rec = NULL, w.don = NULL, type.pred = "random",
             out.pred = FALSE, ...)
}
\arguments{
\item{data.rec}{dataframe including the Xs (predictors, listed in \code{match.vars}) and \code{y.rec} (response; target variable in this dataset)}

\item{data.don}{dataframe including the Xs (predictors, listed in  \code{match.vars}) and \code{z.don} (response; target variable in this dataset)}

\item{match.vars}{vector with the names of the Xs variables to be used as predictors (or set in which select the best predictors with lasso) of respectively \code{y.rec} and \code{z.don}}

\item{y.rec}{character indicating the name of Y target variable in \code{data.rec}. It should be a factor.}

\item{z.don}{character indicating the name of Z target variable in \code{data.don}. It should be a factor.}

\item{pred}{character specifying the method used to obtain predictions of both Y and Z. Available methods include
\code{pred = "multinom"} (default) fits two multinomial models (\pkg{nnet} function \code{\link[nnet]{multinom}}) to get predictions with Y and Z as response variables and \code{match. vars} as predictors; \code{pred = "lasso"} uses the lasso method (\R package \pkg{glmnet}, function \code{\link[glmnet]{cv.glmnet}}) and cross-validation to select a subset of \code{match.vars} that are the best predictors of Y and Z, respectively, and then fits the multinomial models with the selected predictors;
\code{pred = "nb"} uses the Naive Bayes classifier to get predictions of Y and Z respectively (\R package \pkg{naivebayes} function \code{\link[naivebayes]{naive_bayes}}); \code{pred = "rf"} fits randomForest to get predictions of both Y and Z (function \code{\link[randomForest]{randomForest}} in \pkg{randomForest}).
}

\item{w.rec}{name of the variable with the weights of the units in \code{data.rec}, if available (default is NULL); the weights, if available, are only used for estimating bounds, not for fitting models.}

\item{w.don}{name of the variable with the weights of the units in \code{data.don}, if available (default is NULL); the weights, if available, are only used for estimating bounds, not for fitting models.}

\item{type.pred}{string specifying how to obtain the predictions of Y and Z. By default, the fitted models return conditional probabilities (or scores), then if \code{type.pred = "random"} (default), the predicted class of Y (Z) is obtained by a random draw with selection probabilities equal to the estimated conditional probabilities (scores); on the contrary, if \code{type.pred = "mostvoted"}, the predicted class is the one with the highest estimated conditional probability (score).}

\item{out.pred}{ Logical. If TRUE (default is FALSE) returns the input datasets with the estimated conditional probabilities (depending on \code{pred} argument), the predicted class for the target variable (Y or Z) in the dataset (depending on \code{type.pred} argument) and the true observed class of Y (or Z).
}
\item{...}{additional arguments, if needed.}

}

\details{
The function evaluates the uncertainty in estimating the contingency table crossing \code{y.rec} (Y) and \code{z.don} (Z) when the two variables are observed in two different samples related to the same target population, but the samples share a number of common predictors. The evaluation of the uncertainty is equivalent to estimating the bounds for each cell in the contingency table where Y and Z intersect; the bounds can be unconditional (Frechet property) or conditional on the predictions of both Y and Z provided by the models fitted according to the \code{pred} argument. This latter way of working avoids many of the drawbacks of obtaining expectations of conditional bounds when conditioning on many X variables, and allows the inclusion of non-categorical predictors. The final estimation of the conditional bounds is provided by the function \code{\link[StatMatch]{Frechet.bounds.cat}}.
}
\value{a list with the following components: 

\code{up.rec} only when \code{out.pred = TRUE} it corresponds to a smaller version of \code{data.rec} with the estimated conditional probabilities for both Y and Z (depending on \code{pred} argument), the predicted class of Y (depending on \code{type.pred} argument),  the predicted class of Z (depending on \code{type.pred} argument), the true observed class of Y and the predictors (argument \code{match.vars}) (and the weights when \code{w.rec} is specified).

\code{up.don} only when \code{out.pred = TRUE} it corresponds to a smaller version of \code{data.don} with the estimated conditional probabilities for both Y and Z (depending on \code{pred} argument), the predicted class of Y (depending on \code{type.pred} argument), the predicted class of Z (depending on \code{type.pred} argument), the true observed class of Z and the predictors (argument \code{match.vars}) (and the weights when \code{w.don} is specified).

\code{p.xx.ini} the estimated relative frequencies in the table crossing predictions of Y and Z; it is estimated after pooling the samples (weighted average of estimates obtained on the separates samples);

\code{p.xy.ini} the estimated table crossing Y and the predictions of both Y and Z estimated from \code{data.rec} (weights are used if provided with the \code{w.rec} argument);

\code{p.xz.ini} the estimated table crossing Z and the predictions of both Y and Z estimated from \code{data.don} (weights are used if provided with the \code{w.don} argument); 

\code{accuracy} the estimated accuracy in predicting respectively Y and Z with the chosen method (argument \code{pred}) and the available predictors (argument \code{match.vars});

\code{bounds} a data.frame whose columns reports the estimated unconditional and conditional bounds for each cell in the contingency table crossing \code{y.rec}(Y) and \code{z.don} (Z);

\code{uncertainty} the uncertainty associated to input data, measured in terms of average width of uncertainty bounds with and without conditioning on the predictions (for further details see \code{\link[StatMatch]{Frechet.bounds.cat}}.
}
\references{
D'Orazio, M., (2024). \emph{Is Statistical Matching feasible?} Note, \url{https://www.researchgate.net/publication/387699016_Is_statistical_matching_feasible}.
}

\author{
Marcello D'Orazio \email{mdo.statmatch@gmail.com} 
}

\seealso{
  \code{\link[StatMatch]{Frechet.bounds.cat}}.
}

\examples{
data(quine, package="MASS") #loads quine from MASS
str(quine)

# split quine in two subsets
set.seed(223344)
lab.A <- sample(nrow(quine), 70, replace=TRUE)
quine.A <- quine[lab.A, 1:3]
quine.B <- quine[-lab.A, 2:4]

# multinomial model and predictions with most-voted criterion
fbp <- Fbounds.pred(data.rec = quine.A, data.don = quine.B, 
                    match.vars = c("Sex", "Age"), 
                    y.rec = "Eth", z.don = "Lrn", 
                    pred = "multinom", type.pred = "mostvoted")

fbp$p.xx.ini # estimated cross-tab of predictions
fbp$bounds # estimated conditional and unconditional bounds
fbp$uncertainty  # estimated uncertainty about Y*Z

# multinomial model and predictions with randomized criterion
fbp <- Fbounds.pred(data.rec = quine.A, data.don = quine.B, 
                    match.vars = c("Sex", "Age"), 
                    y.rec = "Eth", z.don = "Lrn", 
                    pred = "multinom", type.pred = "random")

fbp$p.xx.ini # estimated cross-tab of predictions
fbp$bounds # estimated conditional and unconditional bounds
fbp$uncertainty  # estimated uncertainty about Y*Z
}

