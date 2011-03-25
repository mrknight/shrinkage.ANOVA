\name{mvr.shrink}
\alias{mvr.shrink}
\alias{mvr.predict}


\title{Multivariate Shrinkage Regression}
\description{
  \code{mvr.shrink} implements an analytic shrinkage approach 
  for multivariate linear regression that can be applied
  for "small n, large p" data.
  
  \code{mvr.predict} is a function to predict the response given 
  the (estimated) regression coefficients.
}

\usage{
mvr.shrink(x, y, lambda, lambda.var, w, verbose=TRUE)
mvr.predict(coef, x)
}

\arguments{
  \item{x}{a predictor matrix.}
  \item{y}{a response vector or matrix.  If \code{y} is missing then
           each individual column of \code{x} serves as
	    univariate response vector that is regressed in turn 
	    against the remaining predictors.}
  \item{lambda}{the correlation shrinkage intensity (range 0-1).
                If \code{lambda} is not specified (the default) it is estimated
		using an analytic formula from Sch\"afer and Strimmer (2005) 
		- see details below.  
		For \code{lambda=0} the empirical correlations are recovered.}
  \item{lambda.var}{the variance shrinkage intensity (range 0-1). 
                If \code{lambda.var} is not specified (the default) it is estimated
	        using an analytic formula from Opgen-Rhein and Strimmer (2007)
                - see details below.  
	        For \code{lambda.var=0} the empirical variances are recovered.}	 
  \item{w}{optional: weights for each data point - if not specified uniform weights are assumed
            (\code{w = rep(1/n, n)} with \code{n = nrow(x)}).}	 
  \item{verbose}{output some status messages while computing (default: TRUE)}
  \item{coef}{the matrix of regression coefficients, as output by \code{mvr.shrink()}.}
  
}
\details{
   
   Regression coefficients and partial correlations/variances are interrelated,
   see for example Cox and Wermuth (1993). Specifically, the regression coefficients
   are the partial correlations times the square-root of the ratio of the partial
   variances.
   
   The function \code{mvr.shrink}
   first computes shrinkage estimates of partial correlation and variance following
   the method described in Sch\"afer and Strimmer (2005) and Opgen-Rhein and
   Strimmer (2007), and then calculates from these the respective regression coefficients. 
   
   Note that if the response is multivariate, the the type of regression performed is
   block regression (Cox and Wermuth 1993). 
}
\value{
  \code{mvr.shrink} returns a matrix with the estimated regression coefficients.

  \code{mvr.predict} returns a response matrix.
}
\author{
  Rainer Opgen-Rhein 
  and Korbinian Strimmer (\url{http://strimmerlab.org}).
  
}
\references{
  Cox, D. R., and N. Wermuth. 1993. Linear dependencies represented by chain graphs.
  Statistical Science \bold{8}:204-283.
  
  Opgen-Rhein, R., and K. Strimmer. 2007. Accurate ranking of 
  differentially expressed genes by a distribution-free shrinkage 
  approach.    Statist. Appl. Genet. Mol. Biol. \bold{6}:9.
  (\url{http://www.bepress.com/sagmb/vol6/iss1/art9/})

  Sch\"afer, J., and K. Strimmer. 2005.  A shrinkage approach to large-scale
  covariance estimation and implications for functional genomics. 
  Statist. Appl. Genet. Mol. Biol. \bold{4}:32.
  (\url{http://www.bepress.com/sagmb/vol4/iss1/art32/})
   
}
\seealso{\code{\link{cov.shrink}}}

\examples{
# load corpcor library
library("corpcor")


# example data
data(longley) 
x     = as.matrix(longley[ 1:13, -2])
y     = as.matrix(longley[ 1:13,  2, drop=FALSE]) 
xtest = as.matrix(longley[14:16, -2])
ytest = as.matrix(longley[14:16,  2, drop=FALSE]) 


# least squares regression is recovered if
# shrinkage intensities are all set to zero
lm(y ~ x)
coefs.ols = mvr.shrink(x, y, lambda=0, lambda.var=0)
coefs.ols


# shrinkage regression coefficients 
# note that these are quite different!
coefs.shrink = mvr.shrink(x, y)
coefs.shrink


# prediction
mvr.predict(coefs.ols, xtest)
mvr.predict(coefs.shrink, xtest)


# if no response matrix is given, then each predictor 
# is regressed in turn against all others
coefs = mvr.shrink(longley)
coefs


# connection between partial correlations and regression coefficients
# (note that the mvr.shrink() function is specifically constructed so that 
# this relationship holds for the shrinkage estimates)

beta = coefs[,-1]    # shrinkage regression coefficients w/o intercept
pc = pcor.shrink(longley) # shrunken partial correlations
k = 3
l = 4
# partial correlations as computed by pcor.shrink()
pc[k,l] 
# partial correlations obtained from mvr.shrink() regression coefficients 
sqrt(beta[k,l] * beta[l,k]) * sign(beta[l,k]) 

}
\keyword{multivariate}
