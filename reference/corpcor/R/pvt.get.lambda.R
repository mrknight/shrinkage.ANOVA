### pvt.get.lambda.R  (2010-03-10)
###
###    Non-public function for computing the shrinkage intensity
###    
###
### Copyright 2005-2010 Korbinian Strimmer
###
### This file is part of the `corpcor' library for R and related languages.
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


# returns lambda to be used for shrinkage
pvt.get.lambda = function(x, lambda, w, verbose, type=c("correlation", "variance"), target)
{
  type = match.arg(type)
  
  if (type == "correlation")
  {
     kind = "lambda (correlation matrix):"
     func = pvt.corlambda
     
     # note: x needs to be the *scaled* data matrix
  }
  
  if (type == "variance")
  {
     kind = "lambda.var (variance vector):"
     func = pvt.varlambda
     
     # note: x needs to be the *centered* data matrix
  }
  
   
  # if lambda = 0: don't shrink
  # if lambda > 0: shrink with given lambda
  # if lambda < 0: shrink with estimated lambda  (the default)
 
     
  if (lambda < 0)
  { 
    if (verbose)
    {
      cat(paste("Estimating optimal shrinkage intensity", kind))     
    }
    
    # estimate optimal shrinkage intensity 
    # target: correlations/covariances -> 0  
    lambda = func(x, w, target)

    lambda.estimated = TRUE
      
    if (verbose)
    {
      cat(paste(" ", round(lambda, 4), "\n", sep=""))     
    }
  }
  else
  {
    if (lambda > 1) lambda = 1
    lambda.estimated = FALSE
    
    if (verbose)
    {
      cat(paste("Specified shrinkage intensity", kind, round(lambda, 4), "\n"))     
    }    
  }
  
  if (type == "correlation") 
  {
     return (list(lambda=lambda, lambda.estimated=lambda.estimated))
  }


  if (type == "variance") 
  {
     return (list(lambda.var=lambda, lambda.var.estimated=lambda.estimated))
  }
	       
}

# input:  centered data matrix  (standardization is NOT checked)
#         weights of each data point
#         variance target 
#
pvt.varlambda = function(xc, w, target)
{
  # bias correction factors
  w2 = sum(w*w)           # for w=1/n this equals 1/n   where n=dim(xc)[1]
  h1 = 1/(1-w2)           # for w=1/n this equals the usual h1=n/(n-1)
  h1w2 = w2/(1-w2)        # for w=1/n this equals 1/(n-1)

  zz = xc^2
  q1 = colSums( sweep(zz, MARGIN=1, STATS=w, FUN="*") )
  q2 = colSums( sweep(zz^2, MARGIN=1, STATS=w, FUN="*") ) - q1^2   
  numerator = sum( q2 )
  denominator = sum( (q1 - target/h1)^2 )

  if(denominator == 0) 
    lambda = 1
  else
    lambda = min(1, numerator/denominator * h1w2)
  
  return (lambda)
}


# input:  scaled data matrix  (standardization is NOT checked)
#         weights of each data point
#         target (ignored, assumed to be 0) 
#
#
# note: the fast algorithm in this function is due to Miika Ahdesm\"aki
#
pvt.corlambda = function(xs, w, target)
{
  # bias correction factors
  w2 = sum(w*w)           # for w=1/n this equals 1/n   where n=dim(xs)[1]
  h1w2 = w2/(1-w2)        # for w=1/n this equals 1/(n-1)

  sw = sqrt(w)
  Q1.squared = (crossprod(sweep(xs, MARGIN=1, STATS=sw, FUN="*")))^2
  Q2 = crossprod(sweep(xs^2, MARGIN=1, STATS=sw, FUN="*")) - Q1.squared
  denominator = sum(Q1.squared)-sum(diag(Q1.squared)) 
  numerator = sum(Q2)-sum(diag(Q2))

  if(denominator == 0) 
    lambda = 1
  else
    lambda = min(1, numerator/denominator * h1w2)
  
  return (lambda)
}

