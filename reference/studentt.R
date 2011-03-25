### studentt.R  (2008-11-19)
###
###    Student t statistic and related stuff
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


# difference of means  ("fold change")

diffmean.stat = function (X, L)
{
  FUN = diffmean.fun(L=L)
  score = FUN(X)
  
  return( score )
}

diffmean.fun = function (L)
{
    if (missing(L)) stop("Class labels are missing!")
    
    function(X)
    { 
      tmp = centroids(X, L, var.pooled=FALSE, var.groups=FALSE, shrink=FALSE, verbose=TRUE)
      
      # differences between the two groups
      diff = tmp$means[,1]-tmp$means[,2]
      
      return(diff)
    }
}


# student t statistic

studentt.stat = function (X, L)
{
  FUN = studentt.fun(L=L)
  score = FUN(X)
  
  return( score )
}

studentt.fun = function (L)
{
    if (missing(L)) stop("Class labels are missing!")
 
    function(X)
    { 
      tmp = centroids(X, L, var.pooled=TRUE, var.groups=FALSE, shrink=FALSE, verbose=TRUE)
      
      # differences between the two groups
      diff = tmp$means[,1]-tmp$means[,2]
      
      # standard error of diff
      n1 = tmp$samples[1]
      n2 = tmp$samples[2]
      v =  tmp$var.pooled   
      sd = sqrt( (1/n1 + 1/n2)*v )
      
      # t statistic
      t = diff/sd
      
      return(t)
    }
}



