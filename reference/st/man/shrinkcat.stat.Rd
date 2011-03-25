\name{shrinkcat.stat}
\alias{shrinkcat.stat}
\alias{shrinkcat.fun}
\title{Correlation-Adjusted t Score}

\usage{
shrinkcat.stat(X, L, group.thresh = 1, group.method=c("neighborhood", "graphcluster"), verbose=TRUE)
shrinkcat.fun(L, group.thresh = 1, group.method=c("neighborhood", "graphcluster"), verbose=TRUE)
}

\description{
  \code{shrinkcat.stat} and \code{shrinkcat.fun} compute a shrinkage
   estimate of the ``correlation-adjusted t score''
   of Zuber and Strimmer (2009).
}

\arguments{
  \item{X}{data matrix. Note that the \emph{columns} correspond to variables (``genes'')
           and the \emph{rows} to samples.}
  \item{L}{vector with class labels for the two groups. }
  \item{group.thresh}{controls the grouping of features. 
                       If \code{group.thresh} equals
                     1 no grouping occurs (the default).}
  \item{group.method}{determines the way the grouping of features is performed: either a correlation neighborhood is used,
                       or the clusters of the correlation threshold graph.}
  \item{verbose}{print out some (more or less useful) information during computation.}
}
\details{
  The cat (``correlation-adjusted t'') score is the product of the square root of the
  inverse correlation matrix with a vector of t scores. The cat score thus describes the 
  contribution of each individual feature in separating the two groups, 
  after removing the effect of all other features.

   In Zuber and Strimmer (2009)
  it is shown that the cat score is
  a natural criterion to rank features in the presence of correlation.
  If there is no correlation, the cat score reduces to the usual t score
  (hence in this case the estimate from \code{shrinkcat.stat} equals that from \code{\link{shrinkt.stat}}).

  Additionally, the option \code{group.thresh} allows to specify a correlation group
  for a feature.  In this neighborhood determined by absolute empirical correlation
   the grouped cat score is computed, which is a variant of the Hotelling T2 statistic
    summarizing the total contribution of the whole feature set for distinguishing the two groups.  
 If the neighborhood contains only a single feature the grouped cat score is identical to the cat score.
  
}
\value{
  \code{shrinkcat.stat} returns a vector containing a shrinkage estimate of the
 ``cat score''  for each variable/gene.

  The corresponding \code{shrinkcat.fun} functions return a function that
  computes the cat score when applied to a data matrix
  (this is very useful for simulations).
}

\references{
Zuber, V., and K. Strimmer. 2009.  Gene ranking and biomarker discovery under correlation.
See \url{http://arxiv.org/abs/0902.0751} for publication details.
}
\author{
  Verena Zuber 
  and
  Korbinian Strimmer (\url{http://strimmerlab.org}).
}
\seealso{\code{\link{shrinkt.stat}}, \code{\link{cst.stat}}, \code{\link{lait.stat}}.}

\examples{
# load st library 
library("st")

# prostate data set
data(singh2002)
X = singh2002$x
L = singh2002$y

dim(X)      # 102 6033 
length(L)   # 102


# shrinkage cat statistic
\dontrun{
score = shrinkcat.stat(X, L)
idx = order(abs(score), decreasing=TRUE)
idx[1:10]
# 610  364 1720 3647 3375  332 3282 3991 1557  914

# compute q-values and local false discovery rates
library("fdrtool")
fdr.out = fdrtool(as.vector(score))
sum(fdr.out$qval < 0.05)
sum(fdr.out$lfdr < 0.2)
}

# compared with:

# shrinkage t statistic 
score = shrinkt.stat(X, L)
idx = order(abs(score), decreasing=TRUE)
idx[1:10]
# 610 1720 3940  914  364  332 3647 4331  579 1068

# Student t statistic
score = studentt.stat(X, L)
idx = order(abs(score), decreasing=TRUE)
idx[1:10]
# 610 1720  364  332  914 3940 4546 1068  579 4331

# difference of means ("Fold Change")
score = diffmean.stat(X, L)
idx = order(abs(score), decreasing=TRUE)
idx[1:10]
# 735  610  694  298  698  292  739 3940  702  721
}
\keyword{univar}
