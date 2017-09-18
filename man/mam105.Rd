\name{mam105}
\alias{mam105}
\alias{mam105.relltest}
\alias{mam105.mt}
\alias{mam105.ass}
\docType{data}
\title{Mammal Phylogenetic Analysis for 105 Trees}
\description{
  Phylogenetic analysis of six mammal species for 105 trees.
}
\usage{
data(mam105)

mam105.mt

mam105.ass

mam105.relltest
}
\format{
  mam105.mt is a matrix of size 3414 * 105.  The (i,j) element is the
  site-wise log-likelihood value at site-i for tree-j for i=1,...,3414,
  and j=1,...,105.

  mam105.ass is a list of length 130 for association vectors. The
  components are t1, t2, ..., t105 for trees, and e1, e2, ..., e25 for
  edges.

  mam105.relltest is an object of class \code{"relltest"} of length 130.
}
\details{
  This data considers 105 tree topologies of the six mammals without
  any constraints. Except for this point, everything is the same as
  \code{mam15}.
}

\seealso{\code{\link{mam15}}.}
\examples{
data(mam105)

## show the results for trees and edges
mam105.relltest # print stat, shtest, bootstrap probabilities, and AIC
summary(mam105.relltest) # print AU p-values

\dontrun{
## simpler script to create mam105.trees
mam105.mt <- read.mt("mam105.mt")
mam105.ass <- read.ass("mam105.ass")
mam105.trees <- relltest(mam105.mt,nb=100000)
}

\dontrun{
## script to create mam105.relltest
mam105.mt <- read.mt("mam105.mt")
mam105.ass <- read.ass("mam105.ass")
mam105.relltest <- relltest(mam105.mt,nb=100000,ass=mam105.ass)
}

\dontrun{
## Parallel version of the above script (but different in random seed)
mam105.mt <- read.mt("mam105.mt")
mam105.ass <- read.ass("mam105.ass")
library(parallel)
cl <- makeCluster(40)
mam105.relltest <- relltest(mam105.mt,nb=100000,ass=mam105.ass,cluster=cl)
}


}
\keyword{datasets}
