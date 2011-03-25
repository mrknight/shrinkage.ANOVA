\name{regularizedt}
\alias{studentt.stat}
\alias{studentt.fun}
\alias{efront.stat}
\alias{efront.fun}
\alias{sam.stat}
\alias{sam.fun}
\alias{samL1.stat}
\alias{samL1.fun}
\alias{modt.stat}
\alias{modt.fun}


\title{Various (Regularized) t Statistics}

\usage{
studentt.stat(X, L)
studentt.fun(L)
efront.stat(X, L, verbose=TRUE)
efront.fun(L, verbose=TRUE)
sam.stat(X, L)
sam.fun(L)
samL1.stat(X, L, method=c("lowess", "cor"), plot=FALSE, verbose=TRUE)
samL1.fun(L, method=c("lowess", "cor"), plot=FALSE, verbose=TRUE)
modt.stat(X, L)
modt.fun(L)
}

\description{
  These functions provide a simple interface to a variety of (regularized) t statistics
  that are commonly used in the analysis of high-dimensional case-control studies.
}

\arguments{
  \item{X}{data matrix. Note that the \emph{columns} correspond to variables (``genes'')
           and the \emph{rows} to samples.}
  \item{L}{factor containing class labels for the two groups.  }
  \item{method}{determines how the smoothing parameter is estimated (applies only to improved SAM statistic \code{samL1}).}
  \item{plot}{output diagnostic plot (applies only to improved SAM statistic \code{samL1}).}
  \item{verbose}{print out some (more or less useful) information during computation.}
}
\details{
\code{studentt.*} computes the standard equal variance t statistic.

\code{efront.*} computes the t statistic using the 90 \% rule of Efron et al. (2001).

\code{sam.*} computes the SAM t statistic of Tusher et al. (2001).
Note that this requires the additional installation of the ``samr'' package.

\code{samL1.*} computes the improved SAM t statistic of Wu (2005).
Note that part of the code in this function is based on the R code providec
by B. Wu.

\code{modt.*} computes the moderated t statistic of Smyth (2004).
Note that this requires the additional installation of the ``limma'' package.

All the above statistics are compared relative to each other 
and relative to the shrinkage t statistic in Opgen-Rhein and Strimmer (2007).
}
\value{
  The *.stat functions directly return the respective statistic for each variable.
  
  The corresponding *.fun functions return a function that produces the respective
  statistics when applied to a data matrix (this is very useful for simulations).
}

\references{
 Opgen-Rhein, R., and K. Strimmer. 2007. Accurate ranking of 
 differentially expressed genes by a distribution-free shrinkage 
 approach.  
 Statist. Appl. Genet. Mol. Biol. \bold{6}:9.
 (\url{http://www.bepress.com/sagmb/vol6/iss1/art9/})  
}
\author{
  Rainer Opgen-Rhein 
  and
  Korbinian Strimmer (\url{http://strimmerlab.org}).
}
\seealso{\code{\link{diffmean.stat}}, \code{\link{shrinkt.stat}}, \code{\link{shrinkcat.stat}}.}

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

# student t statistic
score = studentt.stat(X, L)
order(abs(score), decreasing=TRUE)[1:10]
# [1] 11068   724  9990 11387 11310  9985  9996 11046    43    50

# compute q-values and local false discovery rates
library("fdrtool")
fdr.out = fdrtool(score) 
sum( fdr.out$qval < 0.05 )
sum( fdr.out$lfdr < 0.2 )
fdr.out$param


# Efron t statistic (90 \% rule)
score = efront.stat(X, L)
order(abs(score), decreasing=TRUE)[1:10]
# [1]  4790 10979 11068  1022    50   724  5762    43 10936  9939

# sam statistic
# (requires "samr" package)
#score = sam.stat(X, L)
#order(abs(score), decreasing=TRUE)[1:10]
#[1]  4790 10979  1022  5762    35   970    50 11068 10905  2693

# improved sam statistic
#score = samL1.stat(X, L)
#order(abs(score), decreasing=TRUE)[1:10]
#[1]  1  2  3  4  5  6  7  8  9 10
# here all scores are zero!

# moderated t statistic
# (requires "limma" package)
#score = modt.stat(X, L)
#order(abs(score), decreasing=TRUE)[1:10]
# [1]  4790 10979  1022  5762    35    50 11068   970 10905    43

# shrinkage t statistic
score = shrinkt.stat(X, L)
order(abs(score), decreasing=TRUE)[1:10]
#[1] 10979 11068    50  1022   724  5762    43  4790 10936  9939
}
\keyword{univar}
