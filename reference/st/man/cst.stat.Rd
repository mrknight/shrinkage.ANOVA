\name{cst.stat}
\alias{cst.stat}
\alias{cst.fun}
\title{Correlation-Shared t-Statistic}

\usage{
cst.stat(X, L, verbose=TRUE)
cst.fun(L, verbose=TRUE)
}

\description{
  \code{shrinkcat.stat} and \code{shrinkcat.fun} compute 
   the ``correlation-shared'' t-statistic of Tibshirani and Wassermann (2006).
 }

\arguments{
  \item{X}{data matrix. Note that the \emph{columns} correspond to variables (``genes'')
           and the \emph{rows} to samples.}
  \item{L}{vector with class labels for the two groups. }
  \item{verbose}{print out some (more or less useful) information during computation.}
}
\details{
  The correlation-shared t-statistic for a gene is computed as 
   the average of  t-scores correlated with that gene. For mathematical
   details see Tibshirani and Wasserman (2006).
}
\value{
  \code{cst.stat} returns a vector containing correlation-shared t-statistic for each variable/gene.

  The corresponding \code{cst.fun} functions return a function that
  computes the correlation-shared t-statistic when applied to a data matrix
  (this is very useful for simulations).
}

\references{
Tibshirani, R., and L. Wasserman. 2006.  Correlation-sharing for detection of differential gene expression.
See \url{http://arxiv.org/abs/math/0608061} for publication details.
}
\author{

  Korbinian Strimmer (\url{http://strimmerlab.org}).
}
\seealso{\code{\link{shrinkcat.stat}}, \code{\link{lait.stat}}.}

\examples{
# load st library 
library("st")

# prostate data set
data(singh2002)
X = singh2002$x
L = singh2002$y

dim(X)      # 102 6033 
length(L)   # 102

# correlation shared t statistic
\dontrun{
score = cst.stat(X, L)
idx = order(abs(score), decreasing=TRUE)
idx[1:10]
# [1]  610 1720  364  332  914 3940 4546 1068  579 4331
}

# compared with:

# Student t statistic
score = studentt.stat(X, L)
idx = order(abs(score), decreasing=TRUE)
idx[1:10]
# [1]  610 1720  364  332  914 3940 4546 1068  579 4331


# for the same example using the shrinkage cat score see shrinkcat.stat() 

}
\keyword{univar}
