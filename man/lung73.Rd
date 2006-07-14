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
  and Shimodaira (2001), and is found in \code{data(lung)} of
  \pkg{pvclust} package. We reanalyze it, again, by the script shown in
  Examples. The result of \code{pvclust} is stored in
  \code{lung73.pvclust}, and model fitting to bootstrap probabilities is
  stored in \code{lung73.sb}. The result of AU p-values is sometimes
  very different from that of pvclust. For example, \code{pvclust} with default
  parameter value gave AU p-value of 0.70 for Edge-67, but the
  \code{sbfit} gives AU p-value (named "k.3") of 0.97 for the same
  edge. Note that the raw bootstrap probability (i.e., the ordinary bootstrap
  probability with scale=1) is 0.04.

  The AU p-values for all nodes are shown by the \code{summary} method,
\preformatted{
> summary(lung73.sb)

Corrected P-values (percent):
    raw            k.1            k.2            k.3            model   aic      
1   100.00 (0.00)  100.00 (0.00)  100.00 (0.00)  100.00 (0.00)  poly.3   -19.82  
2   100.00 (0.00)  100.00 (0.00)  100.00 (0.00)  100.00 (0.00)  poly.3   -13.84  
3   100.00 (0.00)  100.00 (0.00)  100.00 (0.00)  100.00 (0.00)  poly.3   -18.83  
4   100.00 (0.00)  100.00 (0.00)  100.00 (0.00)  100.00 (0.00)  poly.3   -17.00  
5   100.00 (0.00)  100.00 (0.00)  100.00 (0.00)  100.00 (0.00)  poly.3   -18.85  
6    65.51 (0.48)   66.75 (0.17)   78.55 (0.21)   79.23 (0.26)  sing.3    26.29  
7   100.00 (0.00)  100.00 (0.00)  100.00 (0.00)  100.00 (0.00)  poly.3   -17.98  
8   100.00 (0.00)   99.99 (0.00)  100.00 (0.00)  100.00 (0.00)  poly.3   -18.82  
9   100.00 (0.00)  100.00 (0.00)  100.00 (0.00)  100.00 (0.00)  poly.2   -20.14  
10   99.44 (0.07)   99.47 (0.03)   99.60 (0.03)   99.60 (0.03)  poly.3   -16.18  
11   99.97 (0.02)   99.99 (0.00)   99.99 (0.00)   99.99 (0.00)  poly.3   -11.85  
12  100.00 (0.00)  100.00 (0.00)  100.00 (0.00)  100.00 (0.00)  poly.3   -14.31  
13  100.00 (0.00)  100.00 (0.00)  100.00 (0.00)  100.00 (0.00)  poly.3    -6.63  
14   98.77 (0.11)   98.59 (0.07)   98.99 (0.06)   98.98 (0.06)  poly.3    -9.32  
15  100.00 (0.00)  100.00 (0.00)  100.00 (0.00)  100.00 (0.00)  poly.3   -19.18  
16   99.97 (0.02)   99.98 (0.00)   99.99 (0.00)   99.99 (0.00)  poly.3    -0.78  
17   40.65 (0.49)   43.79 (0.18)   60.17 (0.14)   61.30 (0.22)  sing.3   271.44  
18  100.00 (0.00)  100.00 (0.00)  100.00 (0.00)  100.00 (0.00)  poly.2   -19.55  
19   90.51 (0.29)   89.93 (0.18)   95.02 (0.13)   95.02 (0.13)  poly.3     5.15  
20  100.00 (0.00)  100.00 (0.00)  100.00 (0.00)  100.00 (0.00)  poly.2    -1.58  
21  100.00 (0.00)  100.00 (0.00)  100.00 (0.00)  100.00 (0.00)  poly.3   -11.71  
22   92.23 (0.27)   92.84 (0.16)   96.21 (0.12)   96.22 (0.12)  poly.3    27.54  
23   48.16 (0.50)   50.09 (0.24)   84.76 (0.17)   88.85 (0.25)  sing.3   159.26  
24   91.65 (0.28)   91.34 (0.17)   95.81 (0.15)   95.92 (0.15)  sing.3   -12.90  
25   93.55 (0.25)   93.09 (0.16)   96.20 (0.14)   96.23 (0.15)  sing.3    14.63  
26   46.51 (0.50)   49.46 (0.20)   77.70 (0.17)   79.78 (0.25)  sing.3   215.12  
27   29.30 (0.46)   30.59 (0.19)   71.89 (0.12)   75.80 (0.22)  sing.3   373.05  
28   59.79 (0.49)   59.99 (0.16)   83.21 (0.15)   83.31 (0.15)  poly.3    11.51  
29   58.22 (0.49)   56.46 (0.20)   86.03 (0.17)   87.32 (0.22)  sing.3    36.54  
30   88.75 (0.32)   88.49 (0.21)   98.01 (0.10)   98.34 (0.10)  sing.3    -9.12  
31  100.00 (0.00)  100.00 (0.00)  100.00 (0.00)  100.00 (0.00)  poly.3   -17.23  
32   95.98 (0.20)   96.19 (0.13)   99.19 (0.06)   99.30 (0.06)  sing.3    -4.85  
33   86.24 (0.34)   86.10 (0.20)   94.72 (0.13)   94.74 (0.13)  poly.3    99.60  
34   96.82 (0.18)   97.56 (0.08)   98.81 (0.05)   98.81 (0.05)  poly.2   128.95  
35   79.18 (0.41)   75.47 (0.20)   91.10 (0.15)   91.11 (0.16)  poly.3   450.34  
36   93.64 (0.24)   94.22 (0.14)   97.19 (0.10)   97.18 (0.10)  poly.3     5.34  
37   93.17 (0.25)   92.53 (0.18)   99.18 (0.05)   99.29 (0.05)  sing.3    18.65  
38   52.76 (0.50)   46.10 (0.10)   64.16 (0.11)   64.15 (0.11)  poly.3   834.31  
39   78.30 (0.41)   77.16 (0.21)   93.09 (0.14)   93.12 (0.14)  poly.3    36.74  
40   99.81 (0.04)   99.88 (0.01)  100.00 (0.00)  100.00 (0.00)  sing.3    -5.63  
41   57.70 (0.49)   57.16 (0.22)   86.15 (0.18)   89.29 (0.25)  sing.3    26.10  
42   99.77 (0.05)   99.81 (0.01)   99.96 (0.00)   99.96 (0.00)  poly.2   -10.32  
43   46.38 (0.50)   45.19 (0.09)   58.38 (0.10)   58.34 (0.10)  poly.3   106.30  
44   31.68 (0.47)   31.50 (0.23)   75.25 (0.11)   84.86 (0.18)  sing.3    60.08  
45   96.84 (0.17)   97.16 (0.10)   99.48 (0.03)   99.48 (0.03)  poly.3    15.76  
46   21.25 (0.41)   22.42 (0.22)   83.13 (0.12)   91.27 (0.15)  sing.3   238.29  
47   77.28 (0.42)   74.80 (0.20)   91.15 (0.13)   91.15 (0.13)  poly.2    86.30  
48   26.88 (0.44)   24.92 (0.21)   71.86 (0.18)   85.77 (0.12)  sing.3   704.07  
49   44.42 (0.50)   44.73 (0.15)   64.70 (0.12)   65.41 (0.17)  sing.3    38.42  
50   16.25 (0.37)   17.35 (0.20)   85.64 (0.14)   93.61 (0.12)  sing.3   106.03  
51   94.94 (0.22)   95.40 (0.13)   98.57 (0.07)   98.56 (0.07)  poly.3   219.59  
52   38.94 (0.49)   39.35 (0.12)   70.35 (0.12)   70.50 (0.13)  poly.3   121.80  
53   37.90 (0.49)   34.95 (0.22)   76.71 (0.14)   78.19 (0.24)  sing.3   955.03  
54   54.89 (0.50)   59.17 (0.24)   88.65 (0.20)   91.44 (0.26)  sing.3   333.49  
55   99.37 (0.08)   99.50 (0.03)   99.94 (0.01)   99.94 (0.01)  poly.3    -8.75  
56   62.34 (0.48)   61.85 (0.13)   74.87 (0.12)   74.87 (0.12)  poly.2    39.08  
57   68.00 (0.47)   67.32 (0.22)   93.42 (0.16)   94.20 (0.18)  sing.3    22.04  
58   21.87 (0.41)   19.32 (0.21)   72.78 (0.17)   86.43 (0.14)  sing.3  2123.18  
59   43.40 (0.50)   46.46 (0.21)   74.14 (0.15)   76.74 (0.24)  sing.3   109.80  
60   21.06 (0.41)   17.03 (0.20)   78.57 (0.15)   89.73 (0.14)  sing.3  1579.74  
61   58.78 (0.49)   47.42 (0.16)   61.77 (0.16)   60.13 (0.20)  poly.3  3159.79  
62   96.02 (0.20)   96.39 (0.12)   99.10 (0.05)   99.10 (0.05)  poly.3    32.10  
63   57.77 (0.49)   55.53 (0.15)   81.15 (0.15)   81.24 (0.15)  poly.3   301.46  
64   16.07 (0.37)   12.18 (0.16)   80.38 (0.20)   93.07 (0.10)  sing.3  1081.85  
65    2.82 (0.17)    1.96 (0.05)   74.20 (0.42)   96.88 (0.11)  sing.3   304.24  
66   15.35 (0.36)   15.45 (0.17)   78.48 (0.14)   88.68 (0.16)  sing.3    79.41  
67    3.65 (0.19)    2.41 (0.06)   76.40 (0.41)   96.90 (0.10)  sing.3   288.59  
68   26.55 (0.44)   28.34 (0.16)   78.31 (0.15)   78.67 (0.16)  poly.3   649.54  
69   30.47 (0.46)   29.55 (0.21)   75.36 (0.14)   76.63 (0.25)  sing.3    -9.11  
70   29.22 (0.45)   30.72 (0.21)   74.87 (0.14)   76.63 (0.24)  sing.3   263.28  
71   25.78 (0.44)   26.89 (0.17)   79.63 (0.15)   80.03 (0.16)  poly.3   767.25  
72  100.00 (0.00)  100.00 (0.00)  100.00 (0.00)  100.00 (0.00)  poly.1   -23.55  
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
  98, 13784-13789; dataset is available from
  \url{http://genome-www.stanford.edu/lung_cancer/adeno/}.
}
\references{
Suzuki, R. and Shimodaira, H. (2006).
pvclust: An R package for hierarchical clustering with p-values,
\emph{Bioinformatics}, 22, 1540-1542; available from
\url{http://www.is.titech.ac.jp/~shimo/prog/pvclust/}.
}
\seealso{\code{\link{sbpvclust}}, \code{\link{sbfit.pvclust}}}

\examples{
\dontrun{
## script to create lung73.pvclust and lung73.sb
## multiscale bootstrap resampling of hierarchical clustering
library(pvclust)
data(lung)
sa <- 10^seq(-2,2,length=13) # wider range of scales than pvclust default
lung73.pvclust <- pvclust(lung,r=1/sa,nboot=10000,weight=T) 
lung73.sb <- sbfit(lung73.pvclust) # model fitting
}

\dontrun{
## Parallel version of the above script
## parPvclust took 80 mins using 40 cpu's
library(snow)
library(pvclust)
data(lung)
cl <- makeCluster(40) # launch 40 cpu's
sa <- 10^seq(-2,2,length=13) # wider range of scales than pvclust default
lung73.pvclust <- parPvclust(cl,lung,r=1/sa,nboot=10000,weight=T) 
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
sblegend("topright",z=plot(lung73.sb[[67]])) # plot curve fitting
summary(lung73.sb[[67]]) # print au p-values

}
\keyword{datasets}
