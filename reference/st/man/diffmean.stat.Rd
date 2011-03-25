\name{diffmean.stat}
\alias{diffmean.stat}
\alias{diffmean.fun}
\alias{rankprod.stat}
\alias{rankprod.fun}


\title{Difference of Means (``Fold Change'') and Rank Products Statistic}

\usage{
diffmean.stat(X, L)
diffmean.fun(L)
rankprod.stat(X, L)
rankprod.fun(L)
}

\description{
  These function compute the difference of group means (``fold change'') and the
  related rank products statistic of Breitling et al. (2004).
}

\arguments{
  \item{X}{data matrix. Note that the \emph{columns} correspond to variables (``genes'')
           and the \emph{rows} to samples.}
  \item{L}{factor containing class labels for the two groups.  }
}
\details{
\code{diffmean.*} computes the difference of means (i.e. the fold-change for log-transformed data).

\code{rankprod.*} computes the two-sided rank products statistic, i.e. the geometric mean of the
ranks of the pairwise absolute mean differences (Breitling et al. 2004).  Note that for consistency with the other functions in this package the \emph{complement} of the averaged ranks is returned 
 (i.e. rank 1 becomes \code{ncol(X)}, rank 2 becomes \code{ncol(X)-1}, etc.). 
}
\value{
  The *.stat functions directly return the respective statistic for each variable.
  
  The corresponding *.fun functions return a function that produces the respective
  statistics when applied to a data matrix (this is very useful for simulations).
}

\references{
 Breitling, R., et al. 2004. Rank products: a simple, yet powerful, new method
to detect differentially regulated genes in replicated microarray experiments.
 FEBS Letters \bold{573}:83-9.
}
\author{
  Korbinian Strimmer (\url{http://strimmerlab.org}).

  This function is in part based on code from Henry Wirth.
}
\seealso{\code{\link{studentt.stat}},\code{\link{shrinkt.stat}}.}

\examples{
# load st library 
library("st")

# load Choe et al. (2005) data
data(choedata)
X <- choe2.mat
dim(X) # 6 11475  
L <- choe2.L
L

# L may also contain some real labels
L = c("group 1", "group 1", "group 1", "group 2", "group 2", "group 2")


# difference of means resp. fold change statistic
score = diffmean.stat(X, L)
order(abs(score), decreasing=TRUE)[1:10]
# [1]  4790  6620  1022 10979   970    35  2693  5762  5885     2

# two-sided rank products statistic
score = rankprod.stat(X, L)
order(score, decreasing=TRUE)[1:10]
# [1]  4790  1022 10979  6620    35  2693   970  5762  5885     2

}
\keyword{univar}
