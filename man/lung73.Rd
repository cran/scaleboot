\name{lung73}
\alias{lung73}
\alias{lung73.pvclust}
\alias{lung73.sb}
\docType{data}
\title{Clustering of 73 Lung Tumors}
\description{
  Bootstrapping hierarchical clustering of the DNA microarray data set
  of 73 lung tissue samples each containing 916 observed genes.
}
\usage{
data(lung73)

lung73.pvclust

lung73.sb
}
\format{
  \code{lung73.pvclust} is an object of class \code{"pvclust"}
  defined in \pkg{pvclust} of Suzuki and Shimodaira (2006).

  \code{lung73.sb} is an object of class \code{"scalebootv"} of length
  72.
}
\details{
  The microarray dataset of Garber et al. (2001) is reanalyzed in Suzuki
  and Shimodaira (2006), and is found in \code{data(lung)} of
  the \pkg{pvclust} package. We reanalyze it, again, by the script shown in
  Examples. The result of \code{pvclust} is stored in
  \code{lung73.pvclust}, and model fitting to bootstrap probabilities
  by the \pkg{scaleboot} package
  is stored in \code{lung73.sb}.
  The AU p-values obtained by using the \pkg{scaleboot} package
  are sometimes very different from those obtained by the \pkg{pvclust}
  package. For example, \code{pvclust} with default parameter value gave
  AU p-value of 0.70 for Edge-67, but the
  \code{sbfit} gives AU p-value (named "k.3") of 0.95 for the same
  edge. Note that the raw bootstrap probability (i.e., the ordinary bootstrap
  probability with scale=1) is 0.04.

  The AU p-values for all nodes are shown by the \code{summary} method,
\preformatted{
> summary(lung73.sb[60:70])

Corrected P-values (percent):
   raw          k.1          k.2          k.3          model  aic    
60 20.21 (0.40) 20.29 (0.18) 71.40 (0.20) 78.98 (0.44) sing.3  80.46 
61 58.45 (0.49) 55.08 (0.17) 63.15 (0.24) 56.34 (0.38) poly.3 575.85 
62 95.68 (0.20) 95.92 (0.10) 98.64 (0.10) 98.61 (0.12) poly.3 -12.01 
63 58.31 (0.49) 57.30 (0.17) 82.09 (0.20) 81.74 (0.28) poly.3  20.74 
64 15.81 (0.36) 15.58 (0.16) 75.36 (0.21) 84.86 (0.37) sing.3  71.47 
65  2.96 (0.17)  2.80 (0.07) 76.73 (0.51) 94.88 (0.20) sing.3  33.34 
66 15.75 (0.36) 15.92 (0.16) 78.02 (0.20) 87.98 (0.29) sing.3   7.30 
67  3.63 (0.19)  3.31 (0.07) 77.02 (0.47) 95.10 (0.17) sing.3  25.11 
68 26.20 (0.44) 27.06 (0.17) 83.06 (0.18) 84.90 (0.27) poly.3   8.67 
69 29.49 (0.46) 29.65 (0.17) 75.37 (0.22) 75.83 (0.34) poly.3 -14.09 
70 28.31 (0.45) 29.04 (0.19) 76.62 (0.17) 81.54 (0.37) sing.3   0.99 
}  

  Shown above are four types of p-values as well as selected model and AIC
  values.  "raw" is
  the ordinary bootstrap probability, "k.1" is equivalent to "raw" but
  calculated from the multiscale bootstrap, "k.2" is equivalent to the
  third-order AU p-value of CONSEL, and finally "k.3" is an improved
  version of AU p-value. By default, we use "k.3" when copying back the
  p-values to an object of class \code{"pvclust"}.

  See Examples below for details.
}
\note{
 The microarray
  dataset is not included in \code{data(lung73)}, but it is found in
  \code{data(lung)} of the \pkg{pvclust} package.
}
\source{
  Garber, M. E. et al. (2001)
  Diversity of gene expression in adenocarcinoma of the lung,
  \emph{Proceedings of the National Academy of Sciences},
  98, 13784-13789 (dataset is available from
  \url{http://genome-www.stanford.edu/lung_cancer/adeno/}).
}
\references{
Suzuki, R. and Shimodaira, H. (2006).
pvclust: An R package for hierarchical clustering with p-values,
\emph{Bioinformatics}, 22, 1540-1542 (software is available from
CRAN or
\url{http://stat.sys.i.kyoto-u.ac.jp/prog/pvclust/}).
}
\seealso{\code{\link{sbpvclust}}, \code{\link{sbfit.pvclust}}}

\examples{
\dontrun{
## script to create lung73.pvclust and lung73.sb
## multiscale bootstrap resampling of hierarchical clustering
library(pvclust)
data(lung)
sa <- 9^seq(-1,1,length=13) # wider range of scales than pvclust default
lung73.pvclust <- pvclust(lung,r=1/sa,nboot=10000) 
lung73.sb <- sbfit(lung73.pvclust) # model fitting
}

\dontrun{
## Parallel version of the above script
## parPvclust took 80 mins using 40 cpu's
library(parallel)
library(pvclust)
data(lung)
cl <- makeCluster(40) # launch 40 cpu's
sa <- 9^seq(-1,1,length=13) # wider range of scales than pvclust default
lung73.pvclust <- parPvclust(cl,lung,r=1/sa,nboot=10000) 
lung73.sb <- sbfit(lung73.pvclust,cluster=cl) # model fitting
}

## replace au/bp entries in pvclust object
data(lung73)
lung73.new <- sbpvclust(lung73.pvclust,lung73.sb) # au <- k.3

library(pvclust)
plot(lung73.new) # draw dendrogram with the new au/bp values
pvrect(lung73.new)

## diagnose edges 61,...,69
lung73.sb[61:69] # print fitting details
plot(lung73.sb[61:69]) # plot curve fitting
summary(lung73.sb[61:69]) # print au p-values
## diagnose edge 67
lung73.sb[[67]] # print fitting
plot(lung73.sb[[67]],legend="topleft") # plot curve fitting
summary(lung73.sb[[67]]) # print au p-values

}
\keyword{datasets}
