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
###
### INTERNAL FUNCTIONS
###

######################################################################
### INTERNAL: MODEL FITTING


### model fitting 1
## bp : vector of bootstrap probabilities
## nb : vector of number of replicates
## sa : vector of sigma^2's
## psi : function(beta,s)
## inits : matrix of initial beta's
## mag : vector of magnification factor for beta

sbfit1 <- function(bp,nb,sa,psi,inits,mag=1,
                   method=NULL,control=NULL) {
  ss <- sqrt(sa)
  lik <- function(par)
    likbinom(pnorm(-sapply(sa,function(s) psi(mag*par,s))/ss),bp,nb)
  chkcoef <- function(par) {
    y <- psi(mag*par,check=TRUE)
    if(!is.null(y)) y$par <- y$beta/mag
    y
  }
  fit <- optims(inits,lik,method=method,control=control,chkcoef=chkcoef)
  fit$mag <- mag
  fit
}

### wls fitting 1
## mat : function(beta,sa,mag=1) for design matrix of z-value
## init : initial beta, but only length and name are used

sbwlsfit1 <- function(bp,nb,sa,mat,init,mag=1,tol=1e-10) {
  bz <- -qnorm(bp)
  uu <- !is.infinite(bz) # use only these elements for wls
  if(sum(uu)<1) return(NULL) # no data
  bz <- bz[uu]; bp <- bp[uu]; nb <- nb[uu]; sa <- sa[uu]
  X <- mat(init,sa,mag)
  if(sum(uu)<ncol(X)) return(NULL)  # too few data
  vv <- (1-bp)*bp/(dnorm(bz)^2*nb) # var(bz)
  fit <- lsfit(X,bz,1/vv,intercept=F) # WLS(weighted least squares)
  list(par=fit$coef)
}

## likbinom : minus log-likelihood of binomial distribution
##
## Arguments:
##  pr : a vector of probability parameters
##  bp : a vector of observed probabilities
##  nb : a vector of sample sizes
##
## Value:
##  likbinom returns the minums of the log-likelihood value.
likbinom <- function(pr,bp,nb) {
  bp2 <- 1-bp; pr2 <- 1-pr;
  -sum(nb*(bp*logx(pr)+bp2*logx(pr2)))
}


######################################################################
### INTERNAL: P-VALUE


### corrected p-value 1
## fit : output of optims (includes var, mag)
## psi : function(beta,s,k)
## k : degree for corrected p-value (default: k=1)
## s : sigma^2 for corrected p-value (default: s=1)
## sp : sigma^2 for prediction (default: sp=-1)

sbpv1 <- function(fit,psi,k=1,s=1,sp=-1) {
  pval <- function(par) pnorm(-psi(fit$mag*par,s,k=k,sp=sp))
  pv <- pval(fit$par)
  h <- nderiv(pval,fit$par)
  pe <- sqrtx(h %*% fit$var %*% h)
  list(pv=pv,pe=pe)
}


######################################################################
### INTERNAL: CONFIDENCE INTERVAL

### find a root of y(x) == y0
## assume y(x) is monotone
##
## x0 : a starting value for search
## y0,x,y : to define the equation
## w : weight for lsfit
## tol : tolerance for chekcing monotonicity

sbfindroot <- function(x0,y0,x,y,w=rep(1,length(x)),tol=0) {
  ## determine if increasing or decreasing
  f <- lsfit(x,y,w)
  if(f$coef[2]<0) {
    y <- -y  # modify as increasing
    y0 <- -y0
  }

  ## prepare
  o <- order(x)
  x <- x[o]; y <- y[o]; w <- w[o]
  n <- length(x)
  i0 <- sum(x<=x0)
  if(i0==0 || i0==n) return(NULL)
  i1 <- i0; i2 <- i0+1
  ymin <- y[i1]; ymax <- y[i2]
  
  ## find a range which includes y[i1] <= y0 <= y[i2]
  for(k in 1:1000) { # should terminate, but for safty
    if(y[i1]>y0) {
      i1 <- i1 - 1; if(i1<1) return(NULL)
      if(y[i1] < ymin) ymin <- y[i1]
      else if(y[i1]-ymin>tol) return(NULL)
    } else {
      i1 <- i1 + 1; if(i1>n) return(NULL)
      if(y[i1]>y0) i1 <- i1 - 1
    }
    if(y[i2]<y0) {
      i2 <- i2 + 1; if(i2>n) return(NULL)
      if(y[i2] > ymax) ymax <- y[i2]
      else if(ymax-y[i2]>tol) return(NULL)
    } else {
      i2 <- i2 - 1; if(i2<1) return(NULL)
      if(y[i2]<y0) i2 <- i2 +1
    }
    if(i1>=i2) return(NULL)
    if(y[i1]<=y0 && y0<=y[i2]) break
  }

  ## linear interpolation
  xv <- as.numeric(x[i1] + (y0-y[i1])* (x[i2]-x[i1])/(y[i2]-y[i1]))
  ye <- min(abs(y0-y[c(i1,i2)]))  # absolute error in y
  xe <- min(abs(xv-x[c(i1,i2)]))  # absolute error in x
  re <- ye*sqrt((w[i1]+w[i2])/2)  # relative error in y
  return(list(xv=xv,xe=xe,ye=ye,re=re))
}


######################################################################
### MISC: NUMERICAL FUNCTIONS

## optims: Optimization with Multiple Initial Values
##
## Description:
##   Apply optim to muliple initial values and
##   return the minimum of the values.
##
## Arguments:
##  coefs : n x m matrix. m sets of initial values of  n dimensions.
##  fn : function to be minimized.
##  ... : passed to optim.
##
## Value:
##  a list of the same structure as the optim function.
##  value : the minimum returned value.
##  init : the initial value for the minimum returned value.
optims <- function(coefs,fn,chkcoef=NULL,eps=1e-3,...) {
  ## objective function with constraints
  fn1 <- function(coef) {
    coef0[pos] <- coef
    fn(coef0)
  }
  ## optim1: one-dimensional optimization
  optim1 <- function(coef) {
    fit <- optimize(fn1,c(coef-10,coef+10))
    x0 <- fit$minimum
    f0 <- fn1(x0);
    fp <- fn1(x0+eps); fm <- fn1(x0-eps)
    hess <- ((fp-f0)-(f0-fm))/(eps*eps)
    list(par=x0,value=f0,hessian=hess)
  }
  ## optim2: switch optim1 and optim
  optim2 <- function(coef) {
    if(length(coef)!=length(pos)) stop("internal error")
    if(length(coef)>=2) optim(coef,fn1,hessian=T,...)
    else optim1(coef)
  }

  ## use each column as a initial parameter vector
  coefs <- as.matrix(coefs)
  m <- ncol(coefs); len <- nrow(coefs)
  coef0 <- rep(0,len); pos <- seq(length=len)
  fit0 <- list(par=coefs[,1],value=fn(coefs[,1])+1,init=coefs[,1])
  for(i in 1:m) {
    fit <- optim2(coefs[,i])
    fit$init <- coefs[,i]
    if(fit$value<fit0$value) fit0 <- fit
  }
  ## modify the optimal parameter if necessary
  if(!is.null(chkcoef) && !is.null(y <- chkcoef(fit0$par))) {
    coef0 <- y$par; pos <- which(y$mask)
    ## length(pos)>0 is assumed here...
    fit0 <- optim2(y$par[y$mask])
    coef0[pos] <- fit0$par
    fit0$par <- coef0
    fit0$mask <- y$mask
    fit0$var <- matrix(0,len,len)
    fit0$var[pos,pos] <- solvex(fit0$hessian)
  } else {
    fit0$mask <- rep(T,len)
    fit0$var <- solvex(fit0$hessian)
  }
  names(fit0$par) <- rownames(coefs)

  fit0
}

## solvex: a stable inverse of symmetric matrix
solvex <- function(x,tol=1e-16) {
  e <- eigen(x,sym=T) # eigen values and vectors
  val <- e$val # eigen values
  ##  u <- abs(val) < tol # almost zero values
  u <- val < tol
  ## val[u] <- 0 # set to zero (genaralized inverse)
  val[u] <- NA # set to NA
  ##  val[u] <- 1/tol # set to a large number
  val[!u] <- 1/val[!u] # take the inverse
  a <- e$vec %*% diag(val,nrow=length(val)) %*% t(e$vec)
  dimnames(a) <- dimnames(x)
  a
}

## logx : allow negative x (linear extrapolation for x smaller than 1e-10)
## print(log(1e-10),digits=15)
## [1] -23.0258509299405
logx <- function(x) {
  a <- x<1e-10
  if(sum(a)>0) {
    x[a] <-  x[a]*1e10 - 24.0258509299405
    x[!a] <- log(x[!a])
  } else {
    x <- log(x)
  }
  x
}

## sqrtx : The square root of x allowing negative x values
sqrtx <- function(x) {
  x <- drop(x); sign(x)*sqrt(abs(x)) 
}

## expsx and logsx (exp for large x, linear for small x)
expsx <- function(x) sign(x)*(exp(abs(x))-1)
logsx <- function(x) sign(x)*log(abs(x)+1)

## nderiv: numerical derivative
nderiv <- function(fn,x,eps=1e-4) {
  apply(diag(length(x))*eps,1,
        function(h) 0.5*(fn(x+h)-fn(x-h))/eps)
}

## sebp: standard error of bootstrap probability
sebp <- function(bp,nb) sqrt(bp*(1-bp)/nb)

######################################################################
### MISC: TEXT FORMAT FUNCTIONS

## myformat
myformat <- function(x,s=NULL,digits=6,h=NULL) {
  ux <- is.na(x); x <- round(x,digits)
  x <- format(x,nsmall=digits)
  x[ux] <- "" # print nothing for NA
  y <- format(x,justify="right")
  if(!is.null(s)) {
    us <- is.na(s); s <- round(s,digits)    
    s <- format(s,nsmall=digits)
    s[us] <- "oo" # infinity sign for NA
    lp <- rep("(",length(s)); rp <- rep(")",length(s))
    ## print nothing for which x is NA
    s[ux] <- ""; lp[ux] <- " "; rp[ux] <- " "; 
    y <- paste(y," ",lp,format(s,justify="right"),rp,sep="")
  }
  if(!is.null(h)) y <- format(c(h,y))
  y
}

## catmat
catmat <- function(x,rn=rownames(x),cn=colnames(x),sep="  ",file="") {
  y <- x
  if(!is.null(rn)) y <- cbind(format(rn),y)
  if(!is.null(cn)) {
    if(!is.null(rn)) cn <- c("",cn)
    y <- rbind(cn,y)
    y <- apply(y,2,function(a) format(a))
  }
  apply(y,1,function(a) cat(a,"\n",sep=sep,file=file,append=T))
  invisible(y)
}

## catpval
## internal for print pval
catpval <- function(pv,pe,digits=NULL) {
  op <- sboptions()
  if(is.null(digits)) digits <- op$digits.pval
  if(op$percent) {
    name <- "percent"
    if(missing(pe)) value <- myformat(100*pv,digits=digits)
    else value <- myformat(100*pv,100*pe,digits=digits)    
  } else {
    name <- "probability"
    if(missing(pe)) value <- myformat(pv,digits=digits+2)
    else value <- myformat(pv,pe,digits=digits+2)
  }
  list(name=name,value=value)
}



######################################################################
### MISC: OTHER FUNCTIONS

### transpose of matrix
## output: t(x) for a matrix
trmat <- function(x,col=T) {
  if(is.matrix(x)) t(x)
  else if(col) as.matrix(x) # column vector
  else t(as.matrix(x)) # row vector
}

### matrix as is
## output: x for a matrix
asmat <- function(x,col=F) {
  if(is.matrix(x)) x
  else if(col) as.matrix(x) # column vector
  else t(as.matrix(x)) # row vector
}

### simplist
## output: simplified vector or matrix of list
simplist <- function(a) sapply(a,function(x) x)

### column-wise apply (with same dimensions)
capply <- function(x,fun,...) {
  y <- apply(x,2,fun,...)
  if(is.vector(y))  y <- t(as.matrix(y))
  dimnames(y) <- dimnames(x)
  y
}

### calc max diff
##
## y[i] := max_{j neq i} x[j] - x[i]
##
maxdif <- function(x) {
  i1 <- which.max(x)  # the largest element
  x <- -x + x[i1]
  x[i1] <- -min(x[-i1])  # the second largest value
  x
}

### calc assmaxdif
##
## y[[i]][j] := max_{k neq ass[[i]]} x[k] - x[ass[[i]][j]]
##
assmaxdif <-  function(x,a) {
  y <- vector("list",length(a))
  names(y) <- names(a)
  for(i in seq(along=a))  y[[i]] <- max(x[ -a[[i]] ]) - x[ a[[i]] ]
  y
}

### sum of row vectors
##
## x = matrix (array of row vectors)
## i = indices (for rows)
##
sumrow <- function(x,i) {
  apply(x[i,,drop=F],2,mean)*nrow(x)
}

### weighted sum of row vectors
##
## x = matrix (array of row vectors)
## w = weight vector (for rows)
##
wsumrow <- function(x,w) {
  apply(w*x,2,sum)*nrow(x)/sum(w)
}

### split numbers
## x: vector of interges
## m: number of splits
##
## output: y = matrix of m*length(x)
## such that apply(y,2,sum) = x
## y[i,] is approximately x/m
splitnum <- function(x,m) {
  n <- length(x)
  xm <- floor(x/m)
  y <- t(matrix(xm,n,m)) # initial distribution
  d <- x - apply(y,2,sum) # left over
  i <- 0
  for(j in 1:n) if(d[j]>0) {
    i <- 1+(seq(from=1+i[length(i)],length=d[j])-1) %% m
    y[i,j] <- y[i,j]+1
  }
  y
}

