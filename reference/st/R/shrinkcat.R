### shrinkcat.R  (2010-05-05)
###
###    Shrinkage Estimation of Correlation-Adjusted t Statistic
###
### Copyright 2008-2010 Verena Zuber and Korbinian Strimmer
###
###
### This file is part of the `st' library for R and related languages.
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


shrinkcat.stat = function (X, L, group.thresh = 1,
  group.method=c("neighborhood", "graphcluster"), verbose=TRUE)
{
  FUN = shrinkcat.fun(L=L, group.thresh=group.thresh,
    group.method=group.method, verbose=verbose)
  score = FUN(X)
  
  return( score )
}


shrinkcat.fun = function (L, group.thresh = 1,
  group.method=c("neighborhood", "graphcluster"), verbose=TRUE)
{
    group.method = match.arg(group.method)

    if (missing(L)) stop("Class labels are missing!")
    if (group.thresh > 1 || group.thresh < 0) stop("group.thresh must be chosen from the interval [0,1]")
  
    function(X)
    {
      p = ncol(X)
      n = nrow(X)   

      tmp = centroids(X, L, var.pooled=TRUE, var.groups=FALSE, 
                         powcor.pooled=TRUE, alpha=-1/2, shrink=TRUE, verbose=verbose)
      n1 = tmp$samples[1]
      n2 = tmp$samples[2]
      
      # differences between the two groups
      diff = tmp$means[,1]-tmp$means[,2]

      # standard error of diff
      n1 = tmp$samples[1]
      n2 = tmp$samples[2]
      v =  tmp$var.pooled   
      sd = sqrt( (1/n1 + 1/n2)*v )
      
          
      # t statistic
      t = diff/sd

      # correlation-adjusted statistic
      if (is.null(dim(tmp$powcor.pooled))) # if there is no correlation
        cat = t
      else
        cat = crossprod(tmp$powcor.pooled, t) # decorrelate t

      cat = as.vector(cat)

      if (group.thresh < 1)
      {
        if (verbose)
        {
          cat("Compute grouped cat scores using empirical correlation threshold", group.thresh, "\n")
          if (group.method == "neighborhood") cat("Grouping method: correlation neighborhood\n")
          if (group.method == "graphcluster") cat("Grouping method: clusters of threshold graph\n")
        }

        cat = pvt.groupcat(X, L, cat, group.thresh, group.method)
        attr(cat, "group.thresh") = group.thresh
        attr(cat, "group.method") = group.method
      }

      attr(cat, "lambda.var") = attr(tmp$var.pooled, "lambda.var")
      attr(cat, "lambda") = attr(tmp$powcor.pooled, "lambda")

      return(cat)
    }
}


# determine grouped cat score
# grouping methods:  "neighborhood" = correlation neighborhood (possibly overlapping groups)
#                    "graphcluster" = cluster in threshold graph (disjoint groups)

pvt.groupcat = function(X, L, cat, group.thresh, group.method)
{
  tmp = centroids(X, L, var.pooled=FALSE, var.groups=FALSE, 
          powcor.pooled=TRUE, alpha=1, shrink=FALSE, verbose=FALSE)
  R = tmp$powcor.pooled 
  p = dim(R)[1]	
  gcat = numeric(p)

  if (group.method == "neighborhood")
  {
    for (i in 1:p)
    {
      w = which( abs(R[i,]) >= group.thresh)
      gcat[i] = sqrt( sum( cat[w]^2 ) ) * sign(cat[i])
    }
  }
  if (group.method == "graphcluster")
  {
    bin.mat = ifelse(abs(R) >= group.thresh, 1, 0)
    require("igraph")
    g = graph.adjacency(bin.mat, mode="undirected")
    mem = clusters(g)$membership
    for (i in 1:p)
    {
      w = which(mem == mem[i])
      gcat[i] = sqrt(sum(cat[w]^2)) * sign(cat[i])
    }
  }

  return (gcat)
}

