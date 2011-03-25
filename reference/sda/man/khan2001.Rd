\name{khan2001}
\alias{khan2001}

\title{Childhood Cancer Study of Khan et al. (2001)}


\usage{
data(khan2001)
}

\description{
Gene expression data (2308 genes for 88 samples) from the
microarray study of Khan et al. (2001).  
}

\details{ This data set contains measurements of the gene
expression of 2308 genes for 88 observations:  29 cases of Ewing 
sarcoma (EWS), 11 cases of Burkitt lymphoma (BL), 18 cases of 
neuroblastoma (NB), 25 cases of rhabdomyosarcoma (RMS), 
and 5 other (non-SRBCT) samples.
}

\format{
  \code{khan2001$x} is a 88 x 2308 matrix containing the expression levels. Note that
   rows correspond to samples, and columns to genes. The row names are the original
   image IDs, and the column names the orginal probe labels.

   \code{khan2001$y} is a factor containing the diagnosis for each sample ("BL", "EWS", "NB", "non-SRBCT", "RMS").
 
   \code{khan2001$descr} provides some annotation for each gene.
}

\source{
  The data are described in Khan et al. (2001) and can be
  obtained from 
  \url{http://cbbp.thep.lu.se/pub/Preprints/01/lu_tp_01_06_supp.html}. Note that the values in 
  \code{khan.data$x} are additionally logarithmized (using natural \code{\link[base]{log}}) for normalization.
}

\references{
 Khan et al. 2001.
 Classification and diagnostic prediction of cancers using gene expression 
 profiling and artificial neural networks. Nature Medicine 7:673--679.
}


\examples{
# load sda library
library("sda")

# load full Khan et al (2001) data set
data(khan2001)
dim(khan2001$x) # 88 2308
hist(khan2001$x)
khan2001$y # 5 levels

# data set containing the SRBCT samples
get.srbct = function()
{
  data(khan2001)
  idx = which( khan2001$y == "non-SRBCT" )
  x = khan2001$x[-idx,]
  y = factor(khan2001$y[-idx])
  descr = khan2001$descr[-idx]

  list(x=x, y=y, descr=descr)
}

srbct = get.srbct()
dim(srbct$x)   # 83 2308
hist(srbct$x)
srbct$y # 4 levels
}


\keyword{datasets}
