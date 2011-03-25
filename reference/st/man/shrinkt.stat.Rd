\name{shrinkt.stat}
\alias{shrinkt.stat}
\alias{shrinkt.fun}
\title{The Shrinkage t Statistic}

\usage{
shrinkt.stat(X, L, var.equal=TRUE, verbose=TRUE)
shrinkt.fun(L, var.equal=TRUE, verbose=TRUE)
}

\description{
  \code{shrinkt.stat} and \code{shrinkt.fun} compute the ``shrinkage t'' statistic
   of Opgen-Rhein and Strimmer (2007).
}

\arguments{
  \item{X}{data matrix. Note that the \emph{columns} correspond to variables (``genes'')
           and the \emph{rows} to samples.}
  \item{L}{factor containing class labels for the two groups. }
  \item{var.equal}{assume equal (default) or unequal variances in each group.}
  \item{verbose}{print out some (more or less useful) information during computation.}
}
\details{
  The ``shrinkage t'' statistic is similar to the usual t statistic, with the
  replacement of the sample variances by corresponding shrinkage estimates.
  These are derived in a distribution-free fashion and with little a priori
  assumptions.   Using the ``shrinkage t'' statistic procduces highly accurate rankings -
  see Opgen-Rhein and Strimmer (2007).  

  The``shrinkage t'' statistic can be generalized to include gene-wise correlation,
  see \code{\link{shrinkcat.stat}}.
}
\value{
  \code{shrinkt.stat} returns a vector containing the ``shrinkage t'' 
  statistic for each variable/gene.

  The corresponding \code{shrinkt.fun} functions return a function that
  produces the ``shrinkage t'' statistics when applied to a data matrix
  (this is very useful for simulations).
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
\seealso{\code{\link{studentt.stat}},
         \code{\link{diffmean.stat}},
         \code{\link{shrinkcat.stat}}.}

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

# shrinkage t statistic (equal variances)
score = shrinkt.stat(X, L)
order(abs(score), decreasing=TRUE)[1:10]

# [1] 10979 11068    50  1022   724  5762    43  4790 10936  9939
#  lambda.var (variance vector):  0.3882


# shrinkage t statistic (unequal variances)
score = shrinkt.stat(X, L, var.equal=FALSE)
order(abs(score), decreasing=TRUE)[1:10]

# [1] 11068    50 10979   724    43  1022  5762 10936  9939  9769
#  lambda.var (variance vector):  0.3673   0.3362


# compute q-values and local false discovery rates
library("fdrtool")
fdr.out = fdrtool(score) 
sum( fdr.out$qval < 0.05 )
sum( fdr.out$lfdr < 0.2 )
fdr.out$param

 
}
\keyword{univar}
