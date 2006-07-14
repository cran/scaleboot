##
##  scaleboot: R package for multiscale bootstrap
##  Copyright (C) 2006 Hidetoshi Shimodaira
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 2 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with this program; if not, write to the Free Software
##  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
##
######################################################################
### MAIN: INTERFACE TO OTHER SOFTWARE

### interface to pvclust for hierarchical clustering
##
## library(pvclust) defines class "pvclust"
## x : pvclust object
## y : sbfits object
## y <- sbfits(x)
sbfit.pvclust <- function(x,...)
  sbfit(as.matrix(x$count)/x$nboot,x$nboot,1/x$r,...)
##
## y <- sbfits(x.old)
## x.new <- mbpvclust(x.old,y) # use "k.2"
## x.new <- mbpvclust(x.old,y,k=3) # use "k.3"
sbpvclust <- function(x,mbs,k=3,k.bp=1,...) {
  a <- summary(mbs,k=c(k.bp,k),...)
  pv.bp <- sapply(a,function(s) s$best$pv[1])
  pe.bp <- sapply(a,function(s) s$best$pe[1])
  pv <- sapply(a,function(s) s$best$pv[2])
  pe <- sapply(a,function(s) s$best$pe[2])
  x$edges[,"au"] <- pv
  x$edges[,"se.au"] <- pe
  x$edges[,"bp"] <- pv.bp
  x$edges[,"se.bp"] <- pe.bp
  x
}

### interface to CONSEL for phylogenetic inference
##
## read "mt" format of consel
##  (a matrix of site-wise log-likelihood values)
##
## ntree: number of trees (i.e., number of items)
## nsite: number of sites (i.e., sample size)
## data: x[i,j] = tree-i, site-j
##
## output: x = matrix of size nsite*ntree
##
## Details of file format:
##
## int ntree, nsite;
## double x[1,1], x[1,2],...,x[1,nsite];
## double x[2,1], x[2,2],...,x[2,nsite];
## ...
## double x[ntree,1], x[ntree,2],...,x[ntree,nsite];
##
read.mt <- function(file,tlab="t") {
  a <- scan(file,comment="#",quiet=T)
  ntree <- a[1]
  nsite <- a[2]
  if(length(a) != (2+ntree*nsite))
    stop(length(a),"!= 2+",ntree,"*",nsite)
  cat("Read",ntree,"items of length",nsite,"\n")
  x <- matrix(a[-c(1,2)],nsite,ntree)
  if(!is.null(tlab)) colnames(x) <- paste(tlab,seq(ncol(x)),sep="")  
  x
}
## read "ass" format of consel
##   (association between trees and edges)
##
## nedge: number of edges (i.e., number of hypotheses)
## ntree: number of trees (i.e., number of items)
## association of i-th edge:
## x[i] = {x[i][1],...,x[i][nt[i]]}
## (x[i][j] is a tree-id starting from 0, but adjust to be from 1)
## association of i-th tree:
## y[i] = {y[i][1],...,y[i][ne[i]]}
## (y[i][j] is a edge-id starting from 0, but adjust to be from 1)
##
## output: list(x,y)
##
## Details of file format:
##
## int nedge;
## int nt[1],x[1][1],...,x[1][nt[1]];
## int nt[2],x[2][1],...,x[2][nt[2]];
## ...
## int nt[nedge],x[nedge][1],...,x[nedge][nt[nedge]];
## int ntree;
## int ne[1],y[1][1],...,y[1][ne[1]];
## int ne[2],y[2][1],...,y[2][ne[2]];
## ...
## int ne[ntree],y[ntree][1],...,y[ntree][ne[ntree]];
##
read.ass <- function(file,identity=TRUE,tlab="t",elab="e") {
  a <- scan(file,comment="#",quiet=T)
  k <- 1
  nedge <- a[k]; k <- k+1
  x <- vector("list",nedge)
  for(i in 1:nedge) {
    nt <- a[k]; k <- k+1
    x[[i]] <- 1+a[seq(k,length=nt)]; k <- k+nt
  }
  if(!is.null(elab)) names(x) <- paste(elab,seq(along=x),sep="")
  ntree <- a[k]; k <- k+1
  y <- vector("list",ntree)
  for(i in 1:ntree) {
    ne <- a[k]; k <- k+1
    y[[i]] <- 1+a[seq(k,length=ne)]; k <- k+ne
  }
  if(!is.null(tlab)) names(y) <- paste(tlab,seq(along=y),sep="")
  if(k-1 != length(a)) stop("size mismatch")
  if(identity) {
    x0 <- vector("list",ntree)
    for(i in 1:ntree) x0[[i]] <- i  # identity associations
    if(!is.null(tlab)) names(x0) <- paste(tlab,seq(along=x0),sep="")
    ans <- c(x0,x) # returns only x without y
  } else ans <- list(x=x,y=y)

  cat("Read",nedge,"items for",ntree,"elements\n")
  ans
}
## read "cnt" format of consel
##   (counts of hypotheses)
##
## ntree: number of trees (or hypotheses, in general)
## nscale: number of scales
## id: tree id's (starting from 0, but adjusted to be from 1)
## lik: observed lik differences of trees
## rs: relative sample sizes (n'/n) ; sa = 1/rs
## nb: numbers of bootstrap replicates
## cnt[i,j]: counts for tree-i, scale-j.
##
## output:
##  list(bps,nb,sa,cnt,id,val)
##   bps = cnt/nb
##   sa = 1/rs
##   val = lik
##
## Details of file format:
##
## int ntree,id[1],id[2],...,id[ntree];
## int ntree; double lik[1],...,lik[ntree];
## int ntree, nscale;
## double rs[1],rs[2],...,rs[nscale];
##  @ repeated ntree times
## int ntree, nscale;
## int nb[1],nb[2],...,nb[nscale];
##  @ repeated ntree times
## int ntree, nscale;
## int cnt[1,1],cnt[1,2],...,cnt[1,nscale];
## int cnt[2,1],cnt[2,2],...,cnt[2,nscale];
## ...
## int cnt[ntree,1],cnt[ntree,2],...,cnt[ntree,nscale];
##
read.cnt <- function(file) {
  a <- scan(file,comment="#",quiet=T)
  k <- 1
  ntree <- a[k]; k <- k+1
  id <- 1+a[seq(k,length=ntree)]; k <- k+ntree
  if(a[k]!=ntree) stop(a[k],"!=",ntree); k <- k+1
  lik <- a[seq(k,length=ntree)]; k <- k+ntree
  if(a[k]!=ntree) stop(a[k],"!=",ntree); k <- k+1
  nscale <- a[k]; k <- k+1
  rs <- a[seq(k,length=nscale)]; k <- k+nscale
  for(i in seq(length=ntree-1)) {
    if(any(a[seq(k,length=nscale)]!=rs)) stop("rs mismatch")
    k <- k+nscale
  }
  if(a[k]!=ntree) stop(a[k],"!=",ntree); k <- k+1
  if(a[k]!=nscale) stop(a[k],"!=",nscale); k <- k+1
  nb <- a[seq(k,length=nscale)]; k <- k+nscale
  for(i in seq(length=ntree-1)) {
    if(any(a[seq(k,length=nscale)]!=nb)) stop("nb mismatch")
    k <- k+nscale
  }
  if(a[k]!=ntree) stop(a[k],"!=",ntree); k <- k+1
  if(a[k]!=nscale) stop(a[k],"!=",nscale); k <- k+1
  cnt <- matrix(a[seq(k,length=nscale*ntree)],nscale,ntree)
  k <- k+nscale*ntree
  if(k-1 != length(a)) stop("size mismatch" ) 
  bps <- cnt/nb
  cnt <- t(cnt); bps <- t(bps)
  
  cat("Read",ntree,"items for",nscale,"scales\n")  
  list(bps=bps,nb=nb,sa=1/rs,cnt=cnt,id=id,val=lik)
}
  
