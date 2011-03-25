\name{lait.stat}
\alias{lait.stat}
\alias{lait.fun}
\alias{lai.tscore}
\title{Correlation-Predicted t-Statistic}

\usage{
lait.stat(X, L, f=0.2, verbose=TRUE)
lait.fun(L, f=0.2, verbose=TRUE)
lai.tscore(gene, tscore, corr, f=0.2, plot=FALSE)
}

\description{
  \code{lait.stat}, \code{laicat.fun}, and \code{lai.tscore}
    compute   the ``correlation-predicted'' t-statistic of Lai (2008).
 }

\arguments{
  \item{X}{data matrix. Note that the \emph{columns} correspond to variables (``genes'')
           and the \emph{rows} to samples.}
  \item{L}{vector with class labels for the two groups. }
  \item{verbose}{print out some (more or less useful) information during computation.}
  \item{f}{smoother span used in \code{\link{lowess}} (default value: 0.2)}
  \item{gene}{the gene for which the Lai t-score is computed}
  \item{tscore}{a vector with t-scores}
  \item{corr}{a matrix containing correlations}
  \item{plot}{show scatter plot correlations versus t-scores with predicted t-score}
  
}
\details{
  The correlation-predicted t-statistic for a gene is 
  the t-score predicted by local linear regression using all other
  genes.  For mathematical details see Lai (2008).
}
\value{
  \code{lait.stat} returns a vector containing correlation-predicted t-statistic for each variable/gene.

  The corresponding \code{lait.fun} functions return a function that
  computes the correlation-shared t-statistic when applied to a data matrix
  (this is very useful for simulations).

  The function \code{lai.tscore} allows to compute the correlation-predicted t-statistic
  for a gene given a correlation matrix and a vector of t-statistics.
}

\references{
Lai, Y.. 2008.  Genome-wide co-expression based prediction of differential expression.
Bioinformatics \bold{24}:666-673.
}
\author{

  Verena Zuber and Korbinian Strimmer (\url{http://strimmerlab.org}).
}
\seealso{\code{\link{shrinkcat.stat}}, \code{\link{cst.stat}}.}

\examples{
# load st library 
library("st")

# prostate data set
data(singh2002)
X = singh2002$x
L = singh2002$y

dim(X)      # 102 6033 
length(L)   # 102

# compute correlation-predicted t-score for various choices
# of smoothing span 

\dontrun{

score1 = lait.stat(X, L, f=0.1)
idx1 = order(abs(score1), decreasing=TRUE)
idx1[1:10]
# 1072  297 1130 4495 4523 4041 1089  955  373 3848

score3 = lait.stat(X, L, f=0.3)
idx3 = order(abs(score3), decreasing=TRUE)
idx3[1:10]
# 1130  962 1688 1223  583 1118  955  297  698 1219

score5 = lait.stat(X, L, f=0.5)
idx5 = order(abs(score5), decreasing=TRUE)
idx5[1:10]
#  698  962 1223 1219  739 1172  583  694 3785 3370 

score7 = lait.stat(X, L, f=0.7)
idx7 = order(abs(score7), decreasing=TRUE)
idx7[1:10]
#  698  739 1219  962 3785  725  694  735 3370 1172


# pick the one with highest correlation to Student t score
t = studentt.stat(X, L)
cor(t, score1, method="spearman") # 0.4265832
cor(t, score3, method="spearman") # 0.471273
cor(t, score5, method="spearman") # 0.4750564
cor(t, score7, method="spearman") # 0.4666669

}

# focus on gene 19
t = studentt.stat(X, L)
R = centroids(X, L,  powcor.pooled=TRUE, alpha=1, 
              shrink=FALSE, verbose=TRUE)$powcor.pooled

lai.tscore(gene=19, t, R, f=0.5, plot=TRUE)
}
\keyword{univar}
