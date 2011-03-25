### shrinkt.R  (2008-10-27)
###
###    Shrinkage t Statistic
###
### Copyright 2006-2008 Rainer Opgen-Rhein and Korbinian Strimmer
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


shrinkt.stat = function (X, L, var.equal=TRUE, verbose=TRUE)
{
  FUN = shrinkt.fun(L=L, var.equal=var.equal, verbose=verbose)
  score = FUN(X)
  
  return( score )
}


shrinkt.fun = function (L, var.equal=TRUE, verbose=TRUE)
{
    if (missing(L)) stop("Class labels are missing!")
  
    function(X)
    {
      p = ncol(X)
      n = nrow(X)   

      if (var.equal) # compute pooled variance
      {
        tmp = centroids(X, L, var.pooled=TRUE, var.groups=FALSE, shrink=TRUE, verbose=verbose)
        n1 = tmp$samples[1]
        n2 = tmp$samples[2]
      
        # differences between the two groups
        diff = tmp$means[,1]-tmp$means[,2]

        # standard error of diff
        n1 = tmp$samples[1]
        n2 = tmp$samples[2]
        v =  tmp$var.pooled   
        sd = sqrt( (1/n1 + 1/n2)*v )
      }
      else # allow different variances in each class
      {
        tmp = centroids(X, L, var.pooled=FALSE, var.groups=TRUE, shrink=TRUE, verbose=verbose)
        n1 = tmp$samples[1]
        n2 = tmp$samples[2]
      
        # differences between the two groups
        diff = tmp$means[,1]-tmp$means[,2]

        v1 = as.vector(tmp$var.groups[,1])
        v2 = as.vector(tmp$var.groups[,2])
   
        # standard error of diff 
        sd = sqrt( v1/n1 + v2/n2 )
      }
          
      # t statistic
      t = diff/sd

      return(t)
    }
}
