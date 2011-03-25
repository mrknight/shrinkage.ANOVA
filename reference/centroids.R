### centroids.R  (2009-02-27)
###
###    Group centroids, variances, and correlations
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



centroids = function(x, L, mean.pooled=FALSE, var.pooled=TRUE, 
  var.groups=FALSE, powcor.pooled=FALSE, alpha=1, shrink=FALSE, verbose=TRUE)
{
  if (!is.matrix(x)) stop("Input x must be a matrix!")
  p = ncol(x)
  n = nrow(x)

  if (length(L) != n) stop("For each sample there must be a class label!")   
  
  groups = pvt.groups(L)
  samples = groups$samples
  cl.count = groups$cl.count
  cl.names = groups$cl.names

  if (verbose)
  {
    cat("Number of variables:", p, "\n")
    cat("Number of observations:", n, "\n")
    cat("Number of classes:", cl.count, "\n\n")
  }

  # means
  mu = array(0, dim=c(p, cl.count))
  colnames(mu) = cl.names
  rownames(mu) = colnames(x)
 
  if(var.pooled || powcor.pooled)
  {
     xc = array(0, dim=c(n,p))  # storage for centered data
     #colnames(xc) = colnames(x)
  }

  if(var.groups)
  {
    # storage for variances
    v = array(0, dim=c(p, cl.count))
    colnames(v) = c(cl.names)
    rownames(v) = colnames(x)
    if (shrink) attr(v, "lambda.var") = numeric(cl.count)
  }
  else
  {
    v = NULL
  }

  for (k in 1:cl.count)
  {
     idx = groups$idx[,k]
     Xk = x[ idx, ,drop = FALSE]
     mu[,k] = colMeans(Xk)

     if(var.pooled || powcor.pooled)
       xc[idx,] = sweep(Xk, 2, mu[,k]) # center data

     if (var.groups)
     {
       if(verbose) cat("Estimating variances (class #", k, ")\n", sep="")
        if (shrink)
        {
          vs = var.shrink(Xk, verbose=verbose)
          v[,k] = as.vector(vs)
          attr(v, "lambda.var")[k] = attr(vs, "lambda.var")
        }
        else
          v[,k] = as.vector(var.shrink(Xk, lambda.var=0, verbose=FALSE))
     }
  }

  if (var.pooled)
  {
    if (verbose) cat("Estimating variances (pooled across classes)\n")
   
    if (shrink)
    {
      v.pool = var.shrink(xc, verbose=verbose)
    }
    else
    {
      v.pool = as.vector(var.shrink(xc, lambda.var=0, verbose=FALSE))
      attr(v.pool, "lambda.var") = NULL
    }
    attr(v.pool, "class") = NULL
    attr(v.pool, "lambda.var.estimated") = NULL
    
    v.pool = v.pool*(n-1)/(n-cl.count)
    names(v.pool) = colnames(x)
  }
  else
  {
    v.pool = NULL
  }
  
  if (powcor.pooled)
  {
    if (verbose)
    {
       if (alpha==1)
         cat("Estimating correlation matrix (pooled across classes)\n")
       else if (alpha==-1)
         cat("Estimating inverse correlation matrix (pooled across classes)\n")
       else
         cat("Estimating correlation matrix to the power of", alpha, "(pooled across classes)\n")
    }

    if (shrink)
    {
      powr = powcor.shrink(xc, alpha=alpha, collapse=TRUE, verbose=verbose)
    }
    else
    {
      powr = powcor.shrink(xc, alpha=alpha, lambda=0, collapse=TRUE, verbose=FALSE)
      attr(powr, "lambda") = NULL
    }
    attr(powr, "class") = NULL
    attr(powr, "lambda.estimated") = NULL      

    # note there is no correction factor for correlation
    if (is.matrix(powr))
    {
      colnames(powr) = colnames(x)
      rownames(powr) = colnames(x)
    }
    else
    {
      names(powr) = colnames(x)
    }
  }
  else
  {
    powr = NULL
  }

  if (mean.pooled)
  {
    mu.pooled = colMeans(x)
  }
  else
  {
    mu.pooled = NULL
  }

  return( list(samples=samples, means=mu, mean.pooled=mu.pooled, 
    var.pooled=v.pool, var.groups=v, powcor.pooled=powr, alpha=alpha))
}



## private function ##

pvt.groups = function(L)
{
   y = factor(L)  # note that this creates new levels (in contrast to as.factor)
     
   cl.names = levels(y)
   cl.count = length(cl.names)

   idx = array(FALSE, dim = c(length(y), cl.count))
   colnames(idx) = cl.names
   nn = integer(cl.count)
   names(nn) = cl.names
   
   for (k in 1:cl.count)
   {
      idx[,k] = ( y == cl.names[k] )
      nn[k] = sum(idx[,k])
   }

   return( list(idx=idx, samples=nn, cl.count=cl.count, cl.names=cl.names) )
}

