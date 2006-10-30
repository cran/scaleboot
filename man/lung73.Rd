\name{lung73}
\alias{lung73}
\alias{lung73.pvclust}
\alias{lung73.sb}
\docType{data}
\title{Clustering of 73 Lung Tumors}
\description{
  Bootstrapping hierarchical clustering of DNA microarray dataset
  of 73 lung tissues with 916 observations of genes.
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
  73.
}
\details{
  The microarray dataset of Garber et al. (2001) is reanalyzed in Suzuki
  and Shimodaira (2006), and is found in \code{data(lung)} of
  \pkg{pvclust} package. We reanalyze it, again, by the script shown in
  Examples. The result of \code{pvclust} is stored in
  \code{lung73.pvclust}, and model fitting to bootstrap probabilities is
  stored in \code{lung73.sb}. The result of AU p-values is sometimes
  very different from that of pvclust. For example, \code{pvclust} with default
  parameter value gave AU p-value of 0.70 for Edge-67, but the
  \code{sbfit} gives AU p-value (named "k.3") of 0.95 for the same
  edge. Note that the raw bootstrap probability (i.e., the ordinary bootstrap
  probability with scale=1) is 0.04.

  The AU p-values for all nodes are shown by the \code{summary} method,
\preformatted{
> summary(lung73.sb)

Corrected P-values (percent):
   raw           k.1           k.2           k.3           model  aic    
1  100.00 (0.00) 100.00 (0.00) 100.00 (0.00) 100.00 (0.00) raw      0.00 
2  100.00 (0.00) 100.00 (0.00) 100.00 (0.00) 100.00 (0.00) poly.3 -18.79 
3  100.00 (0.00) 100.00 (0.00) 100.00 (0.00) 100.00 (0.00) poly.1 -23.62 
4  100.00 (0.00) 100.00 (0.00) 100.00 (0.00) 100.00 (0.00) raw      0.00 
5  100.00 (0.00) 100.00 (0.00)  99.92 (0.57)  99.72 (1.95) poly.3 -18.23 
6   65.99 (0.47)  65.80 (0.20)  79.62 (0.19)  81.93 (0.43) sing.3  -8.37 
7  100.00 (0.00) 100.00 (0.00) 100.00 (0.00) 100.00 (0.00) poly.3 -19.44 
8  100.00 (0.00) 100.00 (0.00) 100.00 (0.00) 100.00 (0.00) poly.2 -20.47 
9  100.00 (0.00) 100.00 (0.00) 100.00 (0.00) 100.00 (0.00) raw      0.00 
10  99.59 (0.06)  99.51 (0.03)  99.63 (0.03)  99.63 (0.03) poly.2 -14.73 
11  99.97 (0.02)  99.98 (0.00)  99.98 (0.01)  99.98 (0.01) poly.3 -16.70 
12 100.00 (0.00) 100.00 (0.00) 100.00 (0.00) 100.00 (0.00) poly.1 -22.78 
13 100.00 (0.00) 100.00 (0.00) 100.00 (0.00) 100.00 (0.00) poly.1 -22.32 
14  98.65 (0.12)  98.59 (0.05)  99.05 (0.06)  99.05 (0.06) poly.2 -11.63 
15 100.00 (0.00) 100.00 (0.00) 100.00 (0.00) 100.00 (0.00) poly.2 -21.23 
16  99.94 (0.02)  99.96 (0.01)  99.97 (0.01)  99.97 (0.01) poly.2 -13.58 
17  41.28 (0.49)  41.41 (0.21)  61.27 (0.19)  65.35 (0.47) sing.3  -5.26 
18 100.00 (0.00) 100.00 (0.00) 100.00 (0.00) 100.00 (0.00) poly.1 -23.26 
19  90.14 (0.30)  90.19 (0.13)  95.42 (0.11)  95.42 (0.11) poly.2 -15.10 
20 100.00 (0.00) 100.00 (0.00) 100.00 (0.00) 100.00 (0.00) poly.2 -18.57 
21 100.00 (0.00) 100.00 (0.00) 100.00 (0.00) 100.00 (0.00) poly.2 -17.25 
22  92.49 (0.26)  92.47 (0.12)  94.93 (0.22)  94.67 (0.28) poly.3   2.53 
23  48.13 (0.50)  48.65 (0.22)  84.94 (0.14)  90.76 (0.25) sing.3 -13.04 
24  91.73 (0.28)  91.26 (0.12)  95.14 (0.20)  95.02 (0.25) poly.3  -0.13 
25  93.34 (0.25)  93.17 (0.11)  96.48 (0.11)  96.48 (0.11) poly.2  -1.34 
26  46.11 (0.50)  46.57 (0.22)  77.72 (0.16)  83.46 (0.35) sing.3  13.25 
27  28.79 (0.45)  29.82 (0.19)  75.69 (0.17)  81.01 (0.37) sing.3   4.12 
28  60.64 (0.49)  60.58 (0.17)  83.03 (0.20)  82.73 (0.28) poly.3  -7.05 
29  59.12 (0.49)  57.67 (0.17)  85.96 (0.18)  86.14 (0.25) poly.3  15.44 
30  89.32 (0.31)  88.67 (0.14)  97.78 (0.11)  97.86 (0.12) poly.3  -4.81 
31 100.00 (0.00) 100.00 (0.00) 100.00 (0.00) 100.00 (0.00) poly.2 -19.07 
32  96.16 (0.19)  96.09 (0.10)  99.05 (0.09)  99.10 (0.12) sing.3  -9.32 
33  86.00 (0.35)  85.73 (0.13)  92.36 (0.22)  91.84 (0.29) poly.3 -11.98 
34  96.87 (0.17)  96.70 (0.09)  96.51 (0.23)  95.99 (0.30) poly.3   2.11 
35  79.04 (0.41)  77.29 (0.15)  87.99 (0.22)  86.35 (0.31) poly.3 251.76 
36  94.16 (0.23)  94.16 (0.11)  97.37 (0.15)  97.47 (0.18) poly.3  -3.68 
37  92.58 (0.26)  92.54 (0.12)  98.98 (0.04)  98.98 (0.04) poly.2  36.09 
38  52.78 (0.50)  52.00 (0.17)  60.86 (0.25)  57.48 (0.37) poly.3  79.06 
39  79.00 (0.41)  77.81 (0.15)  92.15 (0.18)  91.67 (0.24) poly.3  41.84 
40  99.70 (0.05)  99.76 (0.02)  99.98 (0.01)  99.97 (0.01) poly.3  -8.00 
41  56.87 (0.50)  57.83 (0.17)  87.09 (0.17)  88.45 (0.22) poly.3   4.79 
42  99.73 (0.05)  99.76 (0.02)  99.95 (0.01)  99.95 (0.01) poly.2 -11.28 
43  45.38 (0.50)  46.83 (0.17)  58.19 (0.25)  57.61 (0.37) poly.3  -0.82 
44  31.95 (0.47)  31.28 (0.20)  73.25 (0.20)  84.21 (0.28) sing.3  15.25 
45  96.89 (0.17)  96.86 (0.09)  99.22 (0.07)  99.18 (0.08) poly.3  -1.48 
46  21.55 (0.41)  21.29 (0.17)  80.37 (0.19)  92.04 (0.18) sing.3  -3.81 
47  77.47 (0.42)  75.52 (0.15)  89.68 (0.20)  88.93 (0.27) poly.3  56.00 
48  26.97 (0.44)  26.67 (0.19)  67.32 (0.25)  80.74 (0.28) sing.3  59.49 
49  44.39 (0.50)  45.08 (0.17)  66.62 (0.23)  67.35 (0.34) poly.3  -8.08 
50  16.72 (0.37)  16.63 (0.16)  84.23 (0.18)  94.37 (0.15) sing.3 -10.86 
51  94.88 (0.22)  94.43 (0.11)  96.12 (0.20)  95.44 (0.28) poly.3   9.20 
52  38.31 (0.49)  39.70 (0.19)  72.47 (0.20)  73.27 (0.40) sing.3   0.82 
53  37.20 (0.48)  36.36 (0.17)  70.08 (0.23)  68.57 (0.36) poly.3  19.53 
54  54.74 (0.50)  55.86 (0.22)  89.18 (0.13)  95.22 (0.16) sing.3  -1.07 
55  99.55 (0.07)  99.52 (0.03)  99.94 (0.01)  99.94 (0.01) poly.2 -11.44 
56  61.88 (0.49)  62.19 (0.12)  75.11 (0.16)  75.11 (0.16) poly.2   8.15 
57  67.69 (0.47)  67.66 (0.17)  93.64 (0.14)  94.03 (0.17) poly.3  -6.30 
58  21.57 (0.41)  21.43 (0.18)  64.15 (0.25)  78.34 (0.33) sing.3  30.88 
59  43.42 (0.50)  43.74 (0.22)  73.59 (0.17)  80.03 (0.36) sing.3 -14.05 
60  20.21 (0.40)  20.29 (0.18)  71.40 (0.20)  78.98 (0.44) sing.3  80.46 
61  58.45 (0.49)  55.08 (0.17)  63.15 (0.24)  56.34 (0.38) poly.3 575.85 
62  95.68 (0.20)  95.92 (0.10)  98.64 (0.10)  98.61 (0.12) poly.3 -12.01 
63  58.31 (0.49)  57.30 (0.17)  82.09 (0.20)  81.74 (0.28) poly.3  20.74 
64  15.81 (0.36)  15.58 (0.16)  75.36 (0.21)  84.86 (0.37) sing.3  71.47 
65   2.96 (0.17)   2.80 (0.07)  76.73 (0.51)  94.88 (0.20) sing.3  33.34 
66  15.75 (0.36)  15.92 (0.16)  78.02 (0.20)  87.98 (0.29) sing.3   7.30 
67   3.63 (0.19)   3.31 (0.07)  77.02 (0.47)  95.10 (0.17) sing.3  25.11 
68  26.20 (0.44)  27.06 (0.17)  83.06 (0.18)  84.90 (0.27) poly.3   8.67 
69  29.49 (0.46)  29.65 (0.17)  75.37 (0.22)  75.83 (0.34) poly.3 -14.09 
70  28.31 (0.45)  29.04 (0.19)  76.62 (0.17)  81.54 (0.37) sing.3   0.99 
71  25.20 (0.43)  25.95 (0.17)  84.44 (0.18)  85.91 (0.27) poly.3  11.49 
72 100.00 (0.00) 100.00 (0.00) 100.00 (0.00) 100.00 (0.00) raw      0.00 
}  

  There are four types of p-values as well as selected model and AIC
  values are shown above.  "raw" is
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
  \code{data(lung)} of \pkg{pvclust} package.
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
\url{http://www.is.titech.ac.jp/~shimo/prog/pvclust/}).
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
library(snow)
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

\dontrun{
library(pvclust)
plot(lung73.new) # draw dendrogram with the new au/bp values
pvrect(lung73.new)
}

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
