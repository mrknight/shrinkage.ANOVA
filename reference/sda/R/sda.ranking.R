### sda.ranking.R  (2009-03-12)
###
###    Shrinkage discriminant analysis (feature ranking)
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



sda.ranking = function(Xtrain, L, diagonal=FALSE, fdr=TRUE, plot.fdr=FALSE, verbose=TRUE)
{
  if (!is.matrix(Xtrain)) stop("Training data must be given as matrix!")
  if (missing(L)) stop("Class labels are missing!")

  if(verbose) 
  {
    if(diagonal)
      cat("Computing shrinkage t-scores (centroid vs. pooled mean) for feature ranking\n\n")
    else
      cat("Computing shrinkage cat scores (centroid vs. pooled mean) for feature ranking\n\n")
  }

  tmp = centroids(Xtrain, L, mean.pooled=TRUE, var.pooled=TRUE, var.groups=FALSE, 
            powcor.pooled=FALSE, shrink=TRUE, verbose=verbose)
 
  mu = tmp$means
  mup = tmp$mean.pooled
  s2 = tmp$var.pooled
  sc = sqrt(s2)
 
  nk = tmp$samples       # samples in class k
  n = sum(nk)            # number of samples
  p = nrow(mu)           # number of features
  cl.count = ncol(mu)    # number of classes
 
  rm(tmp)


  ############################################################# 
  # compute coefficients for feature ranking
  #############################################################

  # compute cat scores (centroid vs. pooled mean)

  cat = array(0, dim=c(p, cl.count))
  colnames(cat) = paste("cat.", colnames(mu), sep="")
  rownames(cat) = rownames(mu)

  if(!diagonal)
  {
    if(verbose) cat("Computing the square root of the inverse pooled correlation matrix\n")

    invCor.sqrt = centroids(Xtrain, L, mean.pooled=FALSE, var.pooled=FALSE, var.groups=FALSE, 
            powcor.pooled=TRUE, alpha=-1/2, shrink=TRUE, verbose=FALSE)$powcor.pooled

    lambda = attr(invCor.sqrt, "lambda")
    if(verbose) cat("Estimating optimal shrinkage intensity lambda (correlation matrix):", 
                    round(lambda, 4), "\n")
  }

  m = sqrt(1/nk - 1/n) # note the minus sign!
  for (k in 1:cl.count)
  {
    diff = mu[,k]-mup
    if (diagonal)
    { 
      cat[,k] = diff/(m[k]*sc)
    }
    else
    {
      cat[,k] = crossprod(invCor.sqrt, diff/(m[k]*sc) )
    }
  }
  if (!diagonal) rm(invCor.sqrt)

  score = apply(cat^2, 1, sum) # sum of squared cat-scores
  names(score) = rownames(cat)
  idx = order(score, decreasing = TRUE)

  if (fdr)
  {
    if (verbose) cat("\nComputing false discovery rates and higher cricitism scores for each feature\n")

    if (cl.count == 2)
    {
      fdr.out = fdrtool(cat[,1], plot=plot.fdr, verbose=FALSE)
    }
    else
    {
      z = score^(1/3) # Wilson-Hilferty transformation to normality
     
      # center before feeding into fdrtool
      #offset = median(z)
      d = density(z)
      offset = d$x[which.max(d$y)]
      z = z-offset 
      fdr.out = fdrtool(z, plot=plot.fdr, verbose=FALSE)
    }
    lfdr = fdr.out$lfdr # local false discovery rates
    pval = fdr.out$pval # p-values

    # compute absolute HC score for each p-value
    rp = rank(pval)/p
    v = rp*(1-rp)
    v[v==0] = min( v[v > 0] ) # just to make sure we have no zero v
    HC = sqrt(p)*(rp-pval)/sqrt(v)
    HC = abs(HC)

    ranking = cbind(idx, score[idx], cat[idx, , drop=FALSE], lfdr[idx], HC[idx])
    colnames(ranking) = c("idx", "score", colnames(cat), "lfdr", "HC")
    rm(fdr.out)
  }
  else
  {
    ranking = cbind(idx, score[idx], cat[idx, , drop=FALSE])
    colnames(ranking) = c("idx", "score", colnames(cat))
  }
  rm(cat)

  attr(ranking, "class") = "sda.ranking"
  attr(ranking, "diagonal") = diagonal
  attr(ranking, "cl.count") = cl.count

  return(ranking)
}


plot.sda.ranking = function(x, top=40, ...)
{
  if ( class(x) != "sda.ranking" )
    stop ("sda.ranking x needed as input!")

  cl.count = attr(x, "cl.count")
  diagonal = attr(x, "diagonal")
  if (diagonal) 
    xlab = "t-Scores (Centroid vs. Pooled Mean)"
  else 
    xlab = "Correlation-Adjusted t-Scores (Centroid vs. Pooled Mean)"

  top = min( nrow(x), top ) # just to be sure ...

  idx = 2+(1:cl.count)
  cn = colnames(x)[idx]
  colnames(x)[idx] = substr(cn, 5, nchar(cn))

  if (is.null(rownames(x))) rownames(x) = x[, 1]

  score = x[1:top, 2]
  DATA = as.data.frame.table( x[1:top, idx] )

  require("lattice")
  dotplot(reorder(Var1,rep(score, cl.count )) ~ Freq | Var2, 
    data = DATA, origin = 0, type = c("p", "h"), 
    main = paste("The", top, "Top Ranking Features"), 
    xlab = xlab, 
    layout=c(cl.count,1), ...) 
}



