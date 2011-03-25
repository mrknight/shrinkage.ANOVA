### predict.sda.R  (2009-02-27)
###
###    Shrinkage discriminant analysis (prediction)
###
### Copyright 2008-2009 Korbinian Strimmer
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




predict.sda = function(object, Xtest, feature.idx, verbose=TRUE, ...)
{
  if ( missing(object) ) {
    stop("A sda fit object must be supplied.")
  }

  if ( missing(Xtest) ) {
    stop("A new data to predict must be supplied.")
  }
  
  if (!is.matrix(Xtest)) stop("Test data must be given as matrix!")
  ntest = nrow(Xtest)
  freq = object$prior
  cl.count = length(freq)

  if (ncol(Xtest) != nrow(object$predcoef))
    stop("Different number of predictors in sda object (", 
         nrow(object$predcoef), ") and in test data (", 
         ncol(Xtest), ")", sep="")

  if (missing(feature.idx)) feature.idx = 1:nrow(object$predcoef)

  probs = array(0, dim=c(ntest, cl.count) )
  score = numeric(cl.count) 
  yhat = integer(ntest)

  pw = object$predcoef[feature.idx, 1:cl.count+cl.count, drop=FALSE]
  ref = object$predcoef[feature.idx,(1:cl.count), drop=FALSE]

  if (verbose) cat("Prediction uses", nrow(pw), "features.\n")

  for (i in 1:ntest)
  {
    xs = Xtest[i, feature.idx]  # test sample
    for (k in 1:cl.count)
    {
       score[k] = crossprod(pw[, k], xs - ref[, k]) + log(freq[k])
    }
    probs[i,] = score2prob(score)
    yhat[i] = which.max(score)
  }
  probs = zapsmall(probs)
  attr(yhat, "levels") = names(freq)
  class(yhat)= "factor"
  colnames(probs) = names(freq)
  rownames(probs) = rownames(Xtest)

  return(list(class=yhat, posterior=probs) )
}

score2prob = function(x)
{
   x = x-max(x)
   x = exp(x)
   x = x/sum(x)

   return(x)
}

