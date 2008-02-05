##
##  scaleboot: R package for multiscale bootstrap
##  Copyright (C) 2006-2007 Hidetoshi Shimodaira
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
### MAIN: MODEL FITTING

##
## extracting elements
##

"[.scalebootv" <- function(x, i, ...)
  structure(NextMethod("["),
            class = class(x),
            models = attr(x,"models"),
            bps = attr(x,"bps")[i,,drop=F],
            nb = attr(x,"nb"),
            sa = attr(x,"sa")
            )

## general
sbfit <- function(x,...) UseMethod("sbfit")

## refitting
sbfit.scaleboot <- function(x,models=names(x$fi),...) {
  sbfit.default(x$bp,x$nb,x$sa,models,...)
}

## refitting
sbfit.scalebootv <- function(x,models=attr(x,"models"),...) {
  sbfit.matrix(attr(x,"bps"),attr(x,"nb"),attr(x,"sa"),models,...)
}


### default
##
## bp : vector of bootstrap probabilities
## nb : vector of number of replicates
## sa : vector of sigma^2's
## models : model names
##
sbfit.default <- function(x,nb,sa,models=NULL,nofit=FALSE,...) {
  bp <- x
  op <- sboptions()
  if(is.null(models)) models <- op$models
  if(is.numeric(models)) models <- sbmodelnames(m=models)
  bp <- as.vector(bp)
  if(length(bp) != length(sa)) stop("length(bp) != length(sa)")
  nb <- rep(nb,length=length(bp))
  i1 <- order(abs(sa-1))[1] # used for raw bp
  raw <- list(pv=bp[i1],pe=sebp(bp[i1],nb[i1]),nb=nb[i1],s=sa[i1])
  x <- list(bp=bp,nb=nb,sa=sa,raw=raw)
  class(x) <- "scaleboot"
  if(nofit) return(x)
  if(all(bp < 1e-10) || all(bp>1-1e-10)) return(x) # too few data

  ## model fitting
  y <- vector("list",length(models))
  names(y) <- models
  z <- sbname(models)
  lik0 <- likbinom(bp,bp,nb) # lik of the unconstraint model
  for(i in seq(along=models)) {

    ## prepare initial values
    ini <- eval(call(paste("sbini",z$base[[i]],sep="."),
                     z$size[[i]],x,y[seq(1,length=i-1)],z$aux[[i]]))
    ## model fitting
    psi <- paste("sbpsi",z$base[[i]],sep=".")  # psi-function name
    if(op$debug) {
      cat("#################### sbfit\n")
      print(psi)
    }
    ans <- sbfit1(bp,nb,sa,get(psi),
                  ini$inits,ini$mag,ini$omg,ini$trg,
                  method=op$method,control=op$control)
    if(!is.null(ans)) {
      ans$psi <- psi
      
      ## diagnostics of model fitting
      ans$df <- length(bp)-length(ans$par) # degrees of freedom
      ans$rss <- 2*(ans$value-lik0) # likelihood ratio test statistic
      if(ans$df>0) ans$pfit <- pchisq(ans$rss,lower=F,df=ans$df) # p-value
      else ans$pfit <- 1.0
      ans$aic <- ans$rss - 2*ans$df # aic

      ## save
      y[[i]] <- ans
    }
  }
  x$fi <- y
  x
}

## matrix
##
## bps: matrix of bp's (each row is a bp vector)
## nb, sa, models : same as sbfit
## names.hp : hypotheses names
## cluster: for parallel computing (snow package)

sbfit.data.frame <- function(x,...) sbfit(as.matrix(x),...)

sbfit.matrix <- function(x,nb,sa,models=NULL,names.hp=rownames(x),
                         nofit=FALSE,cluster=NULL,...) {
  ## preliminary
  if(is.null(models)) models <- sboptions("models")
  if(is.numeric(models)) models <- sbmodelnames(m=models)
  x <- as.matrix(x)
  nb <- rep(nb,length=ncol(x))
  if(ncol(x) != length(sa)) stop("length(bp) != length(sa")
  nrep <- nrow(x) # number of hypotheses
  if(is.null(names.hp)) names.hp <- seq(length=nrep)
  rownames(x) <- names.hp

  if(nofit) {
    ans <- list()
  } else {
    ## apply sbfit to each bp vector
    arg <- structure(vector("list",nrep),names=names.hp)
    for(i in seq(length=nrep)) arg[[i]] <- x[i,]
    if(is.null(cluster)) ans <- lapply(arg,sbfit,nb,sa,models)
    else ans <- parLapply(cluster,arg,sbfit,nb,sa,models)
    names(ans) <- names.hp
    attr(ans,"models") <- models
  }

  ## class definition
  class(ans) <- "scalebootv"
  attr(ans,"bps") <- x
  attr(ans,"nb") <- nb
  attr(ans,"sa") <- sa
  
  ans
}


###########################
### print

## AIC
sbaic <- function(x,...) UseMethod("sbaic")
sbaic.scaleboot <- function(x,k,...) {
  if(missing(k)) {
    aic <- sapply(x$fi,function(f) f$aic)
  } else {
    aic <- sapply(x$fi,function(f) f$rss - k*f$df)
  }
  if(length(aic)==0) aic <- NULL
  aic
}
sbaic.scalebootv <- function(x,...) {
   lapply(x,sbaic,...)
}


## AIC (left value)
"sbaic<-" <- function(x,value) UseMethod("sbaic<-")
"sbaic<-.scaleboot" <- function(x,value) {
  if(all(names(x$fi)==names(value))) {
    for(i in seq(along=value))
      x$fi[[i]]$aic <- value[[i]]
  } else stop("size mismatch")
  x
}
"sbaic<-.scalebootv" <- function(x,value) {
  for(i in seq(along=value)) {
    sbaic(x[[i]]) <- value[[i]]
  }
  x
}



## extract coefficients of fitting (coef in stats)
coef.scaleboot <- function(object,sd=FALSE,...) {
  fi <- object$fi

  if(is.null(fi)) {
    bv <- be <- NULL
  } else {
    size <- max(sapply(fi,function(f) length(f$par)))
    bv <- be <- matrix(NA,length(fi),size)
    for(i in seq(along=fi)) {
      f <- fi[[i]]
      v <- f$par * f$mag
      e <- sqrt(diag(f$var)) * f$mag
      bv[i,1:length(v)] <- v
      be[i,1:length(e)] <- e
      dimnames(bv) <- dimnames(be) <-
        list(names(fi),paste("beta",0:(size-1),sep=""))
    }
  }
  if(sd) {
    coef <- list(estimate=bv,sd=be)
  } else {
    coef <- bv
  }
  coef
}

coef.scalebootv <- function(object,...) {
  lapply(object,coef,...)
}

## print bp, nb, sa
printbps <- function(bps,nb,sa,digits=NULL) {
  mycatmat <- function(x) catmat(x,cn=seq(ncol(x)),sep=" ")

  if(!is.null(bps)) {
    if(is.matrix(bps)) {
      out <- capply(bps,function(x) catpval(x,digits=digits)$value)
      percent <- catpval(0)$name
    } else {
      a <- catpval(bps)
      percent <- a$name
      out <- matrix(a$value,1)
    }
    cat("\nMultiscale Bootstrap Probabilities (",percent,"):\n",sep="")
    mycatmat(out)
  }

  if(!is.null(nb)) {
    cat("\nNumbers of Bootstrap Replicates:\n")
    mycatmat(matrix(sapply(nb,format),1))
  }

  if(!is.null(sa)) {
    cat("\nScales (Sigma Squared):\n")
    mycatmat(matrix(sapply(sa,format,digits=4),1))
  }
}


## print
print.scaleboot <- function(x,sort.by=c("aic","none"),...) {
  printbps(x$bp,x$nb,x$sa)
  
  fi <- x$fi
  if(is.null(fi)) {
    cat("\nNo Model Fitting\n")
    return(invisible(x))
  }

  ## maximum likelihood estimates of parameters
  a <- coef(x,sd=TRUE)
  betamat <- matrix("",nrow(a$estimate),ncol(a$estimate),
                    dimnames=dimnames(a$estimate))
  for(j in seq(length=ncol(a$estimate))) {
    betamat[,j] <- myformat(c(pi,a$estimate[,j]),c(pi,a$sd[,j]),digits=4)[-1]
  }

  ## goodness of fit
  rss <- sapply(fi,"[[","rss")
  df <- sapply(fi,"[[","df")
  pfit <- sapply(fi,"[[","pfit")
  aic <- sapply(fi,"[[","aic")
  gmat <- cbind(myformat(c(pi,rss),digits=2)[-1],
                format(df),
                myformat(c(pi/10,pfit),digits=4)[-1],
                myformat(c(pi,aic),digits=2)[-1])
  dimnames(gmat) <- list(names(fi),
                         c("rss","df","pfit","aic"))

  ## sort
  sort.by <- match.arg(sort.by)
  j <- switch(sort.by,
              none=1:length(aic),
              aic=order(aic))

  ## print tables
  cat("\nCoefficients:\n")
  catmat(betamat[j,,drop=F])
  cat("\nModel Fitting:\n")
  catmat(gmat[j,,drop=F])

  ## find the best model
  aic0 <- min(aic)
  i <- which(aic==aic0)[1]
  model <- names(fi)[[i]]
  cat("\nBest Model: ",model,"\n")

  invisible(x)
}


print.scalebootv <- function(x,...) {
  ## basic information
  bps <- attr(x,"bps")
  nb <- attr(x,"nb")
  sa <- attr(x,"sa")

  printbps(bps,nb,sa,digits=0)

  ## aic table
  models <- attr(x,"models")
  if(length(x)>0) {
    aics <- matrix(NA,length(x),length(models))
    for(i in seq(along=x)) if(!is.null(x[[i]]$fi)) {
      aic <- sapply(x[[i]]$fi,"[[","aic")
      aics[i,] <- aic[models]
    }
    out <- matrix("",length(x),length(models),
                dimnames=list(names(x),models))
    for(j in seq(along=models)) {
      out[,j] <- format(round(aics[,j],digits=2))
    }
    cat("\nAIC values of Model Fitting:\n")
    catmat(out)
  }

  invisible(x)
}
