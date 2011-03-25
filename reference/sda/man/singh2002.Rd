\name{singh2002}
\alias{singh2002}

\title{Prostate Cancer Study of Singh et al. (2002)}


\usage{
data(singh2002)
}

\description{
Gene expression data (6033 genes for 102 samples) from the
microarray study of Singh et al. (2002).  
}

\details{ This data set contains measurements of the gene
expression of 6033 genes for 102 observations:  52 prostate cancer patients
and 50 healty men.
}

\format{
  \code{singh2002$x} is a 102 x 6033 matrix containing the expression levels.
   The rows contain the samples and the columns the genes.

   \code{singh2002$y} is a factor containing the diagnosis for each sample ("cancer" or "healthy").
}

\source{
  The data are described in Singh et al. (2001) and are provided in exactly the form
  as used by Efron (2008) - see \url{http://www-stat.stanford.edu/~ckirby/brad/papers/Ebaydata.R}. 
}

\references{
 D. Singh et al. 2002. Gene expression correlates of clinical prostate cancer behavior.
Cancer Cell 1:203--209.

Efron, B. 2008. Empirical Bayes estimates for large-scale prediction problems.
Technical Report, Standford University.
}


\examples{
# load sda library
library("sda")

# load Singh et al (2001) data set
data(singh2002)
dim(singh2002$x) # 102 6033
hist(singh2002$x)
singh2002$y # 2 levels
}


\keyword{datasets}
