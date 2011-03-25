### sda.R  (2009-02-27)
###
###    Shrinkage discriminant analysis (training the classifier)
###
### Copyright 2008-09 Miika Ahdesmaki and Korbinian Strimmer
###
###
### This file is part of the `sda' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 3, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA


sda = function(Xtrain, L, diagonal=FALSE, verbose=TRUE)
{
  if (!is.matrix(Xtrain)) stop("Training data must be given as matrix!")
  if (missing(L)) stop("Class labels are missing!")

  # shrinkage intensities
  regularization = rep(NA, 3)
  names(regularization) = c("lambda.freqs", "lambda.var", "lambda")
  regularization[3] = 1 # for diagonal=TRUE


  tmp = centroids(Xtrain, L, mean.pooled=TRUE, var.pooled=TRUE, var.groups=FALSE, 
            powcor.pooled=FALSE, shrink=TRUE, verbose=verbose)
 
  mu = tmp$means
  mup = tmp$mean.pooled
  s2 = tmp$var.pooled
  sc = sqrt(s2)
  regularization[2] = attr(s2, "lambda.var")

  nk = tmp$samples       # samples in class k
  n = sum(nk)            # number of samples
  p = nrow(mu)           # number of features
  cl.count = ncol(mu)    # number of classes
 
  rm(tmp)

  
  ############################################################# 
  # compute coefficients for prediction 
  #############################################################

  # class frequencies
  prior = freqs.shrink( nk, verbose=verbose )
  regularization[1] = attr(prior, "lambda.freqs")
  attr(prior, "lambda.freqs") = NULL

  # reference means
  ref = array(0, dim=c(p, cl.count))
  colnames(ref) = paste("ref.", colnames(mu), sep="")
  rownames(ref) = rownames(mu)
  for (k in 1:cl.count)
  {
    ref[,k] = (mu[,k]+mup)/2
  }

  # prediction weights
  pw = array(0, dim=c(p, cl.count) )
  colnames(pw) = paste("pw.", colnames(mu), sep="")
  rownames(pw) = rownames(mu)

  if(!diagonal)
  {
    if(verbose) cat("\nComputing inverse correlation matrix (pooled across classes)\n")
    ic = centroids(Xtrain, L, mean.pooled=FALSE, var.pooled=FALSE, var.groups=FALSE, 
            powcor.pooled=TRUE, alpha=-1, shrink=TRUE, verbose=FALSE)$powcor.pooled

    # check if estimated correlations are zero
    if ( is.null(dim(ic)) ) 
    {
      regularization[3] = 1
      diagonal = TRUE
    }
    else # compute inverse covariance matrix
    {
      regularization[3] = attr(ic, "lambda")
      invS = sweep(sweep(ic, 1, 1/sc, "*"), 2, 1/sc, "*")
    }
    rm (ic)

    if(verbose) cat("Estimating optimal shrinkage intensity lambda (correlation matrix):", 
                    round(regularization[3], 4), "\n")
  }

  for (k in 1:cl.count)
  {
    diff = mu[,k]-mup
    if (diagonal)
    { 
      pw[,k] = diff/s2
    }
    else
    {
      pw[,k] = crossprod(invS, diff)
    }
  }
  if (!diagonal) rm(invS)



  ############################################################# 

  out = list(regularization=regularization, prior=prior, 
             predcoef=cbind(ref, pw))
  class(out)="sda"

  return (out)
}

