\name{choedata}
\alias{choe2.mat}
\alias{choe2.L}
\alias{ choe2.degenes}
\alias{ choe2.mapping}
\alias{ choe2.probe.name}
\alias{choe2.symbol.name}


\title{A Subset of the Choe et al. (2005) ``Golden Spike'' Experiment}

\description{
These data are expression levels for a subset of the genes
investigated in the Choe et al. (2005) ``Golden Spike''
Affymetrix case-control experiment.  
  
From the original data the 2,535 probe sets for \emph{spike-ins with 
ratio 1:1 were removed}, leaving in total 11,475 genes with 
3 replicates per group, and  1,331 known differentially 
expressed genes.
}

\usage{
data(choedata)
}

\format{
\code{choe2.mat} is a matrix of dimension  6 times 11,475. It contains
the samples in its rows and the genes in its columns.

\code{choe2.L} describes the case control-structure of the experiment,
and  \code{choe2.degenes} indicates the known differentially expressed genes.
\code{choe2.symbol.name},
\code{choe2.probe.name}, and
\code{choe2.mapping} provide additional information on the investigated genes.
}

\source{
\url{http://genomebiology.com/2005/6/2/r16/abstract/}, 
\url{http://www.ccr.buffalo.edu/halfon/spike/}
}


\references{
Choe, S. E., M. Boutros, A. M. Michelson, G. M. Church, and M. ~S. Halfon.
2005. Preferred analysis methods for Affymetrix GeneChips revealed
   by a wholly defined control data set.
\emph{Genome Biology } \bold{6}, R16.
}

\examples{
# load st library
library("st")

# load data set
data(choedata)

# 6 samples, 11,475 genes
dim(choe2.mat)

# two groups (case vs. control
choe2.L

# 1,331 differentially expressed genes
sum(choe2.degenes)

# further information on genes
choe2.symbol.name
choe2.probe.name
choe2.mapping 
}
\keyword{datasets}

