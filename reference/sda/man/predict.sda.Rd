\name{predict.sda}
\alias{predict.sda}

\title{Shrinkage Discriminant Analysis 3: Prediction Step}

\description{
  \code{predict.sda} performs class prediction.
}

\usage{
\method{predict}{sda}(object, Xtest, feature.idx, verbose=TRUE, ...)
}

\arguments{
  \item{object}{An \code{sda} fit object obtained from the function \code{sda}.}
  \item{Xtest}{A matrix containing the test data set. Note that 
                the rows correspond to observations and the columns
                to variables.}
  \item{feature.idx}{A vector indicating which features to employ for prediction (if unspecified all
                     features will be used).}
  \item{verbose}{Report shrinkage intensities (sda) and number of used features (predict.sda).}
   \item{...}{Additional arguments for generic predict.}
}

\value{
   \code{predict.sda} predicts class probabilities for each test sample and returns
   a list with two components:

  \item{class}{a factor with the most probable class assignment for each test sample, and}
   \item{posterior}{a matrix containing the respective class posterior probabilities.}
  
}


\author{
  Miiika Ahdesm\"aki and Korbinian Strimmer (\url{http://strimmerlab.org}).
}
\seealso{\code{\link{sda}}, \code{\link{sda.ranking}}.}

\examples{
# see the examples at the "sda" help page

}
\keyword{multivariate}
