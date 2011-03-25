\name{centroids}
\alias{centroids}

\title{Group Centroids, (Pooled) Variances, and Powers of the Pooled Correlation Matrix}

\usage{
centroids(x, L, mean.pooled=FALSE, var.pooled=TRUE, var.groups=FALSE, 
  powcor.pooled=FALSE, alpha=1, shrink=FALSE, verbose=TRUE)
}
\description{
  \code{centroids} computes group centroids and optionally the pooled mean 
   and pooled variance, the group specific variances, and powers of the pooled correlation matrix.
}

\arguments{
  \item{x}{A matrix  containing the  data set. Note that 
                the rows are sample observations and the columns
                are variables.}
  \item{L}{A factor with the group labels. }
  \item{mean.pooled}{Estimate the pooled mean.}
  \item{var.pooled}{Estimate the pooled variances.}
  \item{var.groups}{Estimate all group-specific variances.}
  \item{powcor.pooled}{Estimate pooled correlation matrix (taken to the power of \code{alpha}).}
  \item{alpha}{exponent for the pooled correlation matrix (default: \code{alpha=1}). }
  \item{shrink}{Use empirical or shrinkage estimator.}
  \item{verbose}{Provide some messages while computing.}
}
\details{
If option \code{shrink=TRUE} then the shrinkage estimators 
 \code{\link[corpcor:cov.shrink]{var.shrink}} from Opgen-Rhein and Strimmer (2007)
and  \code{\link[corpcor:cov.shrink]{cor.shrink}} from Sch\"afer and Strimmer (2005) are used.
}
\value{
  \code{centroids} returns a list
  with the following components:
  \item{samples}{a vector containing the samples sizes in each group,}
   \item{means}{the empirical group means,}
   \item{mean.pooled}{the pooled empirical mean,}
   \item{var.pooled}{a vector containing the pooled variances,}
   \item{var.groups}{a matrix containing the group-specific variances, and}
   \item{powcor.pooled}{a matrix containing the pooled correlation matrix to the power of \code{alpha}
                        (if all correlations are zero a vector containing only the  is returned to save space). }
   \item{alpha}{exponent for the pooled correlation matrix.}
}


\author{
  Korbinian Strimmer (\url{http://strimmerlab.org}).
}
\seealso{
\code{\link[corpcor:cov.shrink]{var.shrink}}, 
\code{\link[corpcor]{powcor.shrink}}.}

\examples{
# load sda library
library("sda")

## prepare data set
data(iris) # good old iris data
X = as.matrix(iris[,1:4])
Y = iris[,5]

## estimate centroids and empirical pooled variances
centroids(X, Y)

## show pooled mean
centroids(X, Y, mean.pooled=TRUE)

## compute group-specific variances
centroids(X, Y, var.groups=TRUE)

## and inverse pooled correlation
centroids(X, Y, var.groups=TRUE, powcor.pooled=TRUE, alpha=-1)

## use shrinkage estimator for variances and correlations
centroids(X, Y, var.groups=TRUE, powcor.pooled=TRUE, alpha=-1, shrink=TRUE)
}
\keyword{multivariate}
