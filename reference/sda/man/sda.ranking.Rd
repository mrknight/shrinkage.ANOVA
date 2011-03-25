\name{sda.ranking}
\alias{sda.ranking}
\alias{plot.sda.ranking}

\title{Shrinkage Discriminant Analysis 1: Feature Ranking}

\description{
  \code{sda.ranking} determines a ranking of features by computing cat scores
   between the group centroids and the pooled mean.

  \code{plot.sda.ranking} provides a graphical visualization of the top ranking features..
}

\usage{
sda.ranking(Xtrain, L, diagonal=FALSE, fdr=TRUE, plot.fdr=FALSE, verbose=TRUE)
\method{plot}{sda.ranking}(x, top=40, ...)
}

\arguments{
  \item{Xtrain}{A matrix  containing the training data set. Note that 
                the rows correspond to observations and the columns
                to variables.}
  \item{L}{A factor with the class labels of the training samples. }
  \item{diagonal}{Chooses between LDA (default, \code{diagonal=FALSE}) and DDA (\code{diagonal=TRUE}).}
  \item{fdr}{compute FDR values and HC scores for each feature.}
  \item{plot.fdr}{Show plot with estimated FDR values.}
  \item{verbose}{Print out some info while computing.}
  \item{x}{An "sda.ranking" object -- this is produced by the sda.ranking() function.}
  \item{top}{The number of top-ranking features shown in the plot (default: 40).}
   \item{...}{Additional arguments for generic plot.}
}
\details{
For each feature and centroid a shrinkage cat scores of the mean versus
the pooled mean is computed.  The overall ranking of a feature is determined
by the sum of the squared cat scores across all centroids.
For the diagonal case (LDA) the cat score reduce to the t-score.  Thus in the
two-class diagonal case the feature are simply ranked according to the
(shrinkage) t-scores.

Calling \code{sda.ranking} should be step 1 in a classification analysis.  Steps 2 and 3 are 
\code{\link{sda}} and \code{\link{predict.sda}}

See Ahdesm\"aki and Strimmer (2009) for details.  For the case of two classes
see Zuber and Strimmer (2009). 
}

\value{

   \code{sda.ranking} returns a matrix with the follwing columns:

  \item{idx}{original feature number}
   \item{score}{sum of the squared cat scores - this determines the overall ranking }
   \item{cat}{for each group and feature the cat score of the centroid versus the pooled mean}

  If \code{fdr=TRUE} then additionally local false discovery rate (FDR) values 
  as well as higher criticism (HC) scores are computed for each feature
  (using \code{\link[fdrtool]{fdrtool}}).

}
\references{
Ahdesm\"aki, A., and K. Strimmer. 2009.  Feature selection in "omics" prediction problems 
using cat scores and false non-discovery rate control.
See \url{http://arxiv.org/abs/0903.2003} for publication details.

Zuber, V., and K. Strimmer. 2009.  Gene ranking and biomarker discovery under correlation.
See \url{http://arxiv.org/abs/0902.0751} for publication details.
}

\author{
  Miiika Ahdesm\"aki and Korbinian Strimmer (\url{http://strimmerlab.org}).
}
\seealso{\code{\link{sda}}, \code{\link{predict.sda}}.}

\examples{
# load sda library
library("sda")

################# 
# training data #
#################

# prostate cancer set
data(singh2002)

# training data
Xtrain = singh2002$x
Ytrain = singh2002$y

######################################### 
# feature ranking (diagonal covariance) #
#########################################

# ranking using t-scores (DDA)
ranking.DDA = sda.ranking(Xtrain, Ytrain, diagonal=TRUE)
ranking.DDA[1:10,]

# plot t-scores for the top 40 genes
plot(ranking.DDA, top=40) 

# number of features with local FDR < 0.8 
# (i.e. features useful for prediction)
sum(ranking.DDA[,"lfdr"] < 0.8)

# number of features with local FDR < 0.2 
# (i.e. significant non-null features)
sum(ranking.DDA[,"lfdr"] < 0.2)

# optimal feature set according to HC score
plot(ranking.DDA[,"HC"], type="l")
which.max( ranking.DDA[1:1000,"HC"] ) 


##################################### 
# feature ranking (full covariance) #
#####################################

# ranking using cat-scores (LDA)
ranking.LDA = sda.ranking(Xtrain, Ytrain, diagonal=FALSE)
ranking.LDA[1:10,]

# plot t-scores for the top 40 genes
plot(ranking.LDA, top=40) 

# number of features with local FDR < 0.8 
# (i.e. features useful for prediction)
sum(ranking.LDA[,"lfdr"] < 0.8)

# number of features with local FDR < 0.2 
# (i.e. significant non-null features)
sum(ranking.LDA[,"lfdr"] < 0.2)

# optimal feature set according to HC score
plot(ranking.LDA[,"HC"], type="l")
which.max( ranking.LDA[1:1000,"HC"] ) 

}
\keyword{multivariate}
