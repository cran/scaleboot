\name{mam15}
\alias{mam15}
\alias{mam15.relltest}
\alias{mam15.mt}
\alias{mam15.ass}
\docType{data}
\title{Mammal Phylogenetic Analysis for 15 Trees}
\description{
  Phylogenetic analysis of six mammal species for 15 trees.
}
\usage{
data(mam15)

mam15.mt

mam15.ass

mam15.relltest
}
\format{
  mam15.mt is a matrix of size 3414 * 15.  The (i,j) element is the
  site-wise log-likelihood value at site-i for tree-j for i=1,...,3414,
  and j=1,...,15.

  mam15.ass is a list of length 25 for association vectors. The
  components are t1, t2, ..., t15 for trees, and e1, e2, ..., e10 for
  edges.

  mam15.relltest is an object of class \code{"relltest"} of length 25.
}
\details{
  An example of phylogenetic analysis of six mammal species:
  Homo sapiens (human), Phoca vitulina (harbor seal), Bos taurus (cow),
  Oryctolagus cuniculus (rabbit), Mus musculus (mouse), Didelphis
  virginiana (opossum).  The data is stored in the file \file{mam15.aa},
  which contains amino acid sequences of length N=3414 for the six
  species obtained from mtDNA (see Note below).
  We consider 15 tree topologies of the six mammals as stored in the
  file \file{mam15.tpl}; 

\preformatted{  
((Homsa,(Phovi,Bosta)),Orycu,(Musmu,Didvi)); t1
(Homsa,Orycu,((Phovi,Bosta),(Musmu,Didvi))); t2
(Homsa,((Phovi,Bosta),Orycu),(Musmu,Didvi)); t3
(Homsa,(Orycu,Musmu),((Phovi,Bosta),Didvi)); t4
((Homsa,(Phovi,Bosta)),(Orycu,Musmu),Didvi); t5
(Homsa,((Phovi,Bosta),(Orycu,Musmu)),Didvi); t6
(Homsa,(((Phovi,Bosta),Orycu),Musmu),Didvi); t7
(((Homsa,(Phovi,Bosta)),Musmu),Orycu,Didvi); t8
(((Homsa,Musmu),(Phovi,Bosta)),Orycu,Didvi); t9
(Homsa,Orycu,(((Phovi,Bosta),Musmu),Didvi)); t10
(Homsa,(((Phovi,Bosta),Musmu),Orycu),Didvi); t11
((Homsa,((Phovi,Bosta),Musmu)),Orycu,Didvi); t12
(Homsa,Orycu,(((Phovi,Bosta),Didvi),Musmu)); t13
((Homsa,Musmu),Orycu,((Phovi,Bosta),Didvi)); t14
((Homsa,Musmu),((Phovi,Bosta),Orycu),Didvi); t15
}  

  The log-likelihood values are calculated by PAML software (Ziheng 1997)
  for phylogenetic inference. The two files \file{mam15.aa} and \file{mam15.tpl}
  are fed into PAML to generate the file \file{mam15.lnf} of site-wise
  log-likelihood values.

  Using CONSEL software (Shimodaira and Hasegawa 2001), we convert
  \file{mam15.lnf} and \file{mam15.tpl} to be suitable for scaleboot
  package. Here, we do not use CONSEL for calculating AU p-values, but
  use it only for file conversion. We type

\preformatted{
seqmt --paml mam15.lnf
treeass --outgroup 6 mam15.tpl > mam15.log
}

  The first line above generates \file{mam15.mt}, which is a simple text
  file containing a matrix of site-wise log-likelihood values. The
  second line above generates \file{mam15.ass} and \file{mam15.log}, which
  contain information regarding which edges are included in a tree.  A
  part of \file{mam15.log} is as follows.

\preformatted{  
# leaves: 6
6
  1 Homsa
  2 Phovi
  3 Bosta
  4 Orycu
  5 Musmu
  6 Didvi

# base edges: 10
10 6
          
    123456
  1 +++---  ; (0.2000)
  2 ++++--  ; (0.2000)
  3 +--+--  ; (0.2000)
  4 -+++--  ; (0.2000)
  5 ---++-  ; (0.2000)
  6 +--++-  ; (0.2000)
  7 -++++-  ; (0.2000)
  8 +++-+-  ; (0.2000)
  9 +---+-  ; (0.2000)
 10 -++-+-  ; (0.2000)
}

  The above defines edges named e1,...e10 (base edges) as clusters of six
  mammal species. For example, e1 = +++--- = (Homsa, Phovi, Bosta).

  The converted files are read by scaleboot package in R:

\preformatted{  
mam15.mt <- read.mt("mam15.mt")
mam15.ass <- read.ass("mam15.ass")
}

  \code{mam15.mt} is a matrix of size 3414 * 6 for the site-wise
  log-likelihood values. For testing trees, we need \code{mam15.mt} but
  \code{mam15.ass}.  \code{mam15.ass} is used for testing edges, and is
  a list of length 25 for association vectors for t1,t2,...,t15, and
  e1,e2,...,e10.  For example, \code{mam15.ass$t1 = 1}, indicating tree
  "t1" is included in tree "t1", and \code{mam15.ass$e1 = c(1, 5, 8)},
  indicating edge "e1" is included in trees "t1", "t5", and "t8".

  Multiscale bootstrap resampling is performed by the function
  \code{relltest}. The simplest way to get AU p-values for trees is:

\preformatted{
mam15.trees <- relltest(mam15.mt) # resampling and fitting
summary(mam15.trees) # calculates AU p-values
}

  The \code{\link{relltest}} returns an object of class \code{"relltest"}.
  It calls the function \code{\link{scaleboot}} internally with
  the number of bootstrap replicates \code{nb=10000}, and takes about 20
  mins. Typically, \code{nb=10000} is large enough, but it would be safe
  to use larger value, say \code{nb=100000} as in the Examples below.
  
  Note that the default value of scales in \code{\link{relltest}} ranges
  much wider than that of CONSEL. It is
  \code{sa=10^seq(-2,2,length=13)} for \code{\link{relltest}}, and
  is \code{sa=1/seq(from=0.5,to=1.4,by=0.1)} for CONSEL.

  The \code{mam15.relltest} object in \code{data(mam15)} is similar
  to \code{mam15.trees} above, but is calculated also for
  edges using \code{mam15.ass}. We can extract the result for trees by
\preformatted{
mam15.trees <- mam15.relltest[1:15]
}
  
  The results for trees stored in \code{mam15.trees} object above are in
  the order specified in the columns of \code{mam15.mt}. To sort it by
  increasing order of log-likelihood difference, we can type
\preformatted{
stat <- attr(mam15.trees,"stat")  # the log-likelihood differences
o <- order(stat) # sort it in increasing order
mam15.trees <- mam15.trees[o] # same as mam15.trees in Examples
}

  The results of fitting are shown by the \code{print} method.
\preformatted{
> mam15.trees

Test Statistic, and Shimodaira-Hasegawa test:
     stat   shtest        
t1   -2.66  94.44 (0.07)  
t3    2.66  80.08 (0.13)  
t2    7.40  57.98 (0.16)  
t5   17.57  17.55 (0.12)  
t6   18.93  14.48 (0.11)  
t7   20.11  11.64 (0.10)  
t4   20.60  11.05 (0.10)  
t15  22.22   7.61 (0.08)  
t8   25.38   3.33 (0.06)  
t14  26.32   3.33 (0.06)  
t13  28.86   1.75 (0.04)  
t9   31.64   0.66 (0.03)  
t11  31.75   0.59 (0.02)  
t10  34.74   0.21 (0.01)  
t12  36.25   0.13 (0.01)  

Multiscale Bootstrap Probabilities (percent):
    1   2  3  4  5  6  7  8  9  10 11 12 13 
t1  100 99 95 87 78 69 58 46 34 27 21 17 13 
t3    0  1  5 13 22 29 32 29 24 19 15 12  9 
t2    0  0  0  0  0  1  4  7 10 12 12 13 14 
t5    0  0  0  0  0  0  1  4  6  7  7  6  5 
t6    0  0  0  0  0  1  3  6  8  9  9  9  9 
t7    0  0  0  0  0  0  1  2  4  6  6  6  6 
t4    0  0  0  0  0  0  1  4  5  6  7  7  6 
t15   0  0  0  0  0  0  0  1  2  3  4  4  5 
t8    0  0  0  0  0  0  0  0  1  2  3  5  5 
t14   0  0  0  0  0  0  0  1  3  4  6  6  7 
t13   0  0  0  0  0  0  0  0  1  2  3  4  4 
t9    0  0  0  0  0  0  0  0  0  1  3  4  6 
t11   0  0  0  0  0  0  0  0  0  1  2  3  5 
t10   0  0  0  0  0  0  0  0  0  0  1  2  3 
t12   0  0  0  0  0  0  0  0  0  0  1  2  3 

Numbers of Bootstrap Replicates:
1     2     3     4     5     6     7     8     9     10    11    12    13    
1e+05 1e+05 1e+05 1e+05 1e+05 1e+05 1e+05 1e+05 1e+05 1e+05 1e+05 1e+05 1e+05 

Scales (Sigma Squared):
1    2       3       4   5      6      7 8     9     10    11    12    13    
0.01 0.02154 0.04642 0.1 0.2154 0.4642 1 2.154 4.639 10.01 21.61 46.14 100.4 

AIC values of Model Fitting:
     poly.1     poly.2    poly.3   sing.3   
t1   219783.17  14953.26  2885.53  2455.02  
t3   228338.40  11285.41  2191.95  3237.38  
t2    98108.84   9314.50  1798.68    10.02  
t5   128773.41   4866.12   842.99    12.76  
t6   129431.99  10889.79  1592.49    53.21  
t7    98052.87   4851.55   467.46    -4.81  
t4   123734.55   7911.70  1553.24    -8.97  
t15   66196.14   3174.49   403.00   -16.31  
t8    30512.50   1330.68   216.40    -9.74  
t14   72262.39   4511.50   794.59   -15.73  
t13   42502.49   1991.00   336.79    -8.00  
t9    25700.97   1970.16   214.58    16.57  
t11   20305.67   1494.20    96.66    28.16  
t10   10091.05    667.74    52.68    16.63  
t12    8677.45    702.91    76.19    47.86  
}

  The AU p-values are shown by the \code{summary} method.
\preformatted{
> summary(mam15.trees)

Corrected P-values (percent):
     raw           k.1           k.2           k.3           model   aic      
t1   57.84 (0.16)  54.37 (0.06)  73.67 (0.05)  75.48 (0.07)  sing.3  2455.02  
t3   31.82 (0.15)  28.91 (0.03)  45.77 (0.05)  45.93 (0.05)  poly.3  2191.95  
t2    3.61 (0.06)   3.57 (0.03)  13.67 (0.22)  18.29 (0.41)  sing.3    10.02  
t5    1.32 (0.04)   1.29 (0.02)   6.53 (0.17)   7.85 (0.24)  sing.3    12.76  
t6    3.11 (0.05)   3.12 (0.03)  14.57 (0.23)  19.34 (0.41)  sing.3    53.21  
t7    0.50 (0.02)   0.52 (0.02)   3.70 (0.16)   4.92 (0.25)  sing.3    -4.81  
t4    1.46 (0.04)   1.44 (0.02)   9.50 (0.25)  13.11 (0.44)  sing.3    -8.97  
t15   0.08 (0.01)   0.07 (0.01)   1.04 (0.10)   1.64 (0.19)  sing.3   -16.31  
t8    0.00 (0.00)   0.00 (0.00)   0.05 (0.02)   0.10 (0.04)  sing.3    -9.74  
t14   0.23 (0.02)   0.22 (0.01)   2.62 (0.19)   4.39 (0.40)  sing.3   -15.73  
t13   0.01 (0.00)   0.01 (0.00)   0.23 (0.05)   0.47 (0.12)  sing.3    -8.00  
t9    0.00 (0.00)   0.00 (0.00)   0.22 (0.02)   1.39 (0.13)  sing.3    16.57  
t11   0.00 (0.00)   0.00 (0.00)   0.07 (0.01)   0.55 (0.07)  sing.3    28.16  
t10   0.00 (0.00)   0.00 (0.00)   0.00 (0.00)   0.01 (0.00)  sing.3    16.63  
t12   0.00 (0.00)   0.00 (0.00)   0.00 (0.00)   0.01 (0.00)  sing.3    47.86  
}

  The p-values for 15 trees are shown above.  "raw" is
  the ordinary bootstrap probability, "k.1" is equivalent to "raw" but
  calculated from the multiscale bootstrap, "k.2" is equivalent to the
  third-order AU p-value of CONSEL, and finally "k.3" is an improved
  version of AU p-value.

  The details for each tree are shown by extracting the element. For
  example, the details for the seventh-largest tree in the
  log-likelihood value ("t4") is obtained by
\preformatted{
mam15.trees[[7]] # same as mam15.trees$t4

Multiscale Bootstrap Probabilities (percent):
1    2    3    4    5    6    7    8    9    10   11   12   13   
0.00 0.00 0.00 0.00 0.01 0.22 1.46 3.56 5.43 6.48 6.75 6.69 5.96 

Numbers of Bootstrap Replicates:
1     2     3     4     5     6     7     8     9     10    11    12    13    
1e+05 1e+05 1e+05 1e+05 1e+05 1e+05 1e+05 1e+05 1e+05 1e+05 1e+05 1e+05 1e+05 

Scales (Sigma Squared):
1    2       3       4   5      6      7 8     9     10    11    12    13    
0.01 0.02154 0.04642 0.1 0.2154 0.4642 1 2.154 4.639 10.01 21.61 46.14 100.4 

Coefficients:
        beta0            beta1            beta2             
poly.1  3.9815 (0.0080)                                     
poly.2  2.4081 (0.0072)  0.1615 (0.0005)                    
poly.3  2.0738 (0.0076)  0.2515 (0.0013)  -0.0012 (0.0000)  
sing.3  1.6697 (0.0129)  0.5175 (0.0082)   0.3058 (0.0072)  

Model Fitting:
        rss        df  pfit    aic        
poly.1  123758.55  12  0.0000  123734.55  
poly.2    7933.70  11  0.0000    7911.70  
poly.3    1573.24  10  0.0000    1553.24  
sing.3      11.03  10  0.3554      -8.97  

Best Model:  sing.3 
> summary(mam15.trees[[7]])

Raw Bootstrap Probability:  1.46 (0.04) 

Corrected P-values (percent):
        k.1          k.2          k.3           aic        
poly.1  0.00 (0.00)  0.00 (0.00)   0.00 (0.00)  123734.55  
poly.2  0.51 (0.01)  1.23 (0.02)   1.23 (0.02)    7911.70  
poly.3  1.01 (0.02)  3.39 (0.06)   3.43 (0.06)    1553.24  
sing.3  1.44 (0.02)  9.50 (0.25)  13.11 (0.44)      -8.97  

Best Model:  sing.3

> sblegend("topleft",z=plot(mam15.trees[[7]])) # plot diagnostics
}  
  Especially, the plot diagnostics in the bottom line is useful
  to identify which model is fitting best.

  See other examples below.
}
\note{
  Dataset files for the phylogenetic inference are found at
  \url{http://www.is.titech.ac.jp/~shimo/prog/scaleboot/}.
  For unix users, download \file{mam15-files.tgz}, and for Windows users
  download \file{mam15-files.zip}.  This dataset
  is originally used in Shimodaira and Hasegawa (1999).

  \file{mam15.aa}: amino acid sequences (N=3414) for the six mammals.
  \file{mam15.ass}: association vectors for edges and trees.
  \file{mam15.cnt}: multiscale bootstrap counts.
  \file{mam15.lnf}: site-wise log-likelihood values (output from PAML).
  \file{mam15.log}: detailed information for the associations.
  \file{mam15.mt}: site-wise log-likelihood values (output from seqmt).
  \file{mam15.tpl}: 15 tree topologies.
  }
\source{
  H. Shimodaira and M. Hasegawa (1999). Multiple comparisons of
  log-likelihoods with applications to phylogenetic inference, \emph{Molecular
  Biology and Evolution}, 16, 1114-1116.
}
\references{
  Yang, Z. (1997). PAML: a program package for phylogenetic analysis by
  maximum likelihood, \emph{Computer Applications in BioSciences},
  13:555-556; available from
  \url{http://abacus.gene.ucl.ac.uk/software/paml.html}.
  
  Shimodaira, H. and Hasegawa, M. (2001). CONSEL: for assessing the
  confidence of phylogenetic tree selection, \emph{Bioinformatics}, 17,
  1246-1247; available from
  \url{http://www.is.titech.ac.jp/~shimo/prog/consel/}.
}
\seealso{\code{\link{relltest}}, \code{\link{summary.scalebootv}},
  \code{\link{read.mt}}, \code{\link{read.ass}}.}
\examples{
data(mam15)

## show the results for trees and edges
mam15.relltest # print stat, shtest, bootstrap probabilities, and AIC
summary(mam15.relltest) # print AU p-values

## extract trees in increasing order of the test statistic
stat <- attr(mam15.relltest,"stat") # log-likelihood difference
o <- order(stat[1:15]) # sort in increasing order (t1,...,t15)
mam15.trees <- mam15.relltest[o] # extract elements from relltest object
mam15.trees  # print stat, shtest, bootstrap probabilities, and AIC
summary(mam15.trees) # print AU p-values
plot(mam15.trees) # draw the curve fitting

## extract edges in increasing order of the test statistics
o <- 15+order(stat[16:25]) # sort in increasing order (e1,...,e10)
mam15.edges <- mam15.relltest[o]
mam15.edges  # print stat, shtest, bootstrap probabilities, and AIC
summary(mam15.edges) # print AU p-values
plot(mam15.edges) # draw the curve fitting

\dontrun{
## simpler script to create mam15.trees
mam15.mt <- read.mt("mam15.mt")
mam15.ass <- read.ass("mam15.ass")
mam15.trees <- relltest(mam15.mt,nb=100000)
}

\dontrun{
## script to create mam15.relltest
## It took 260 mins (Xeon 2.4GHz)
mam15.mt <- read.mt("mam15.mt")
mam15.ass <- read.ass("mam15.ass")
mam15.relltest <- relltest(mam15.mt,nb=100000,ass=mam15.ass)
}

\dontrun{
## Parallel version of the above script (but different in random seed)
## It took 13 mins (40 cpu's of Athlon MP 2000+)
mam15.mt <- read.mt("mam15.mt")
mam15.ass <- read.ass("mam15.ass")
library(snow)
cl <- makeCluster(40)
mam15.relltest <- relltest(mam15.mt,nb=100000,ass=mam15.ass,cluster=cl)
}


}
\keyword{datasets}
