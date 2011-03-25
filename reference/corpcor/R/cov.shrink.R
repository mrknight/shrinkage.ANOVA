### cov.shrink.R  (2008-12-01)
###
###    Shrinkage Estimation of Variance Vector, Correlation Matrix,
###    and Covariance Matrix, and their inverses
###
### Copyright 2005-08 Juliane Schaefer, Rainer Opgen-Rhein and Korbinian Strimmer
###
###
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


# correlation
cor.shrink = function(x, lambda, w, collapse=FALSE, verbose=TRUE)
{
   x = as.matrix(x)
   n = nrow(x)
   if (missing(lambda)) lambda = -1  # estimate correlation shrinkage parameter
   w = pvt.check.w(w, n)
   
   # shrinkage correlation
   r = pvt.powscor(x=x, alpha=1, lambda=lambda, w=w, collapse=collapse, verbose=verbose)
   if (verbose) cat("\n")

   return(r)
}


# inverse correlation
invcor.shrink = function(x, lambda, w, collapse=FALSE, verbose=TRUE)
{
   x = as.matrix(x)
   n = nrow(x)  
   if (missing(lambda)) lambda = -1  # estimate correlation shrinkage parameter
   w = pvt.check.w(w, n)
   
   # inverse shrinkage correlation
   invr = pvt.powscor(x=x, alpha=-1, lambda=lambda, w=w, collapse=collapse, verbose=verbose)
   if (verbose) cat("\n")
   
   return(invr)
}


# variances
var.shrink = function(x, lambda.var, w, verbose=TRUE)
{
  x = as.matrix(x)
  n = nrow(x)  
  if (missing(lambda.var)) lambda.var = -1  # estimate variance shrinkage parameter
  w = pvt.check.w(w, n)
  
  # shrinkage variance 
  sv = pvt.svar(x=x, lambda.var=lambda.var, w=w, verbose=verbose)
  if (verbose) cat("\n")
  
  return(sv)
}


# covariance
cov.shrink = function(x, lambda, lambda.var, w, collapse=FALSE, verbose=TRUE)
{   
   x = as.matrix(x)
   n = nrow(x)   
   if (missing(lambda)) lambda = -1          # estimate correlation shrinkage parameter
   if (missing(lambda.var)) lambda.var = -1  # estimate variance shrinkage parameter  
   w = pvt.check.w(w, n)
   
   # shrinkage scale factors
   sc = sqrt( pvt.svar(x=x, lambda.var=lambda.var, w=w, verbose=verbose) )

   # shrinkage correlation
   c = pvt.powscor(x=x, alpha=1, lambda=lambda, w=w, collapse=collapse, verbose=verbose)
   
   # shrinkage covariance 
   if (is.null(dim(c)))
     c = c*sc*sc
   else
     c = sweep(sweep(c, 1, sc, "*"), 2, sc, "*")

   attr(c, "lambda.var") = attr(sc, "lambda.var")
   attr(c, "lambda.var.estimated") = attr(sc, "lambda.var.estimated")
   if (verbose) cat("\n")
                    
   return(c)
}


# precision matrix (inverse covariance)
invcov.shrink = function(x, lambda, lambda.var, w, collapse=FALSE, verbose=TRUE)
{   
   x = as.matrix(x)
   n = nrow(x) 
   if (missing(lambda)) lambda = -1          # estimate correlation shrinkage parameter
   if (missing(lambda.var)) lambda.var = -1  # estimate variance shrinkage parameter  
   w = pvt.check.w(w, n)

   # shrinkage scale factors
   sc = sqrt( pvt.svar(x=x, lambda.var=lambda.var, w=w, verbose=verbose) )
        
   # inverse shrinkage correlation
   invc = pvt.powscor(x=x, alpha=-1, lambda=lambda, w=w, collapse=collapse, verbose=verbose)
   
   # inverse shrinkage covariance 
   if (is.null(dim(invc)))
     invc = invc/sc/sc
   else
     invc = sweep(sweep(invc, 1, 1/sc, "*"), 2, 1/sc, "*")
   
   attr(invc, "lambda.var") = attr(sc, "lambda.var")
   attr(invc, "lambda.var.estimated") = attr(sc, "lambda.var.estimated")
   if (verbose) cat("\n")
   
   return(invc)
}




