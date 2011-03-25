\encoding{latin1}
\name{sda-package}
\alias{sda-package}
\docType{package}
\title{The sda package}

\description{This package performs linear and diagonal discriminant analysis
with variable selection.

The classifier is trained using James-Stein-type shrinkage estimators.
Variable selection is based on a univariate test statistics
(multi-class cat score) and subsequent false non-discovery rate (FNDR) testing.
In addition, higher criticism (HC) feature scores are computed and can also be used 
for feature thresholding.

This approach is particularly suited for high-dimensional classification.
For details see Ahdesm\"aki and Strimmer (2009).

Typically the functions in this package are applied in three steps:
\itemize{
\item{Step 1:}{feature selection with \code{\link{sda.ranking}},} 
\item{Step 2:}{training the classifier with \code{\link{sda}}, and}
\item{Step 3:}{classification using \code{\link{predict.sda}}.} 
}
The accompanying web site (see below) provides example R scripts to illustrate
the functionality of this package.
}

\author{Miika Ahdesm\"aki and Korbinian Strimmer (\url{http://strimmerlab.org/})}
\references{
  See website: \url{http://strimmerlab.org/software/sda/}
}
\keyword{multivariate}
\seealso{
 \code{\link{sda.ranking}}, \code{\link{sda}}, \code{\link{predict.sda}}. 
}
