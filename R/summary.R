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
### MAIN: SUMMARY


##
## extracting elements
##

"[.summary.scalebootv" <- function(x, i, ...)
  structure(NextMethod("["),
            pvalues = attr(x,"pvalues")
            )
##
## scaleboot
##

summary.scaleboot <- function(object,models=names(object$fi),
                              k=1:3,s=1,sp=-1,...) {
  class(object) <- c("summary.scaleboot",class(object))
  pvnames <- paste("k",k,sep=".")

  ## models
  if(is.numeric(models)) models <- names(object$fi)[models]
  
  ## restrict models
  a <- match(models,names(object$fi))
  a <- a[!is.na(a)]
  if(length(a)>0) object$fi <- object$fi[a] else object$fi <- NULL
  models <- names(object$fi)

  ## no model fitting
  if(is.null(object$fi)) {
    pv <- rep(object$raw$pv,length(pvnames))
    pe <- rep(object$raw$pe,length(pvnames))
    names(pv) <- names(pe) <- pvnames
    object$best <- list(model="raw",aic=0,pv=pv,pe=pe)
    return(object)
  }

  ## corrected p-values
  pv <- pe <- matrix(NA,length(models),length(k),
                     dimnames=list(models,pvnames))
  for(i in seq(along=models)) {
    m <- models[[i]]
    psi <- get(object$fi[[m]]$psi)
    for(j in seq(along=k)) {
      y <- sbpv1(object$fi[[m]],psi,k=k[[j]],s=s,sp=sp)
      pv[i,j] <- y$pv
      pe[i,j] <- y$pe
    }
  }
  object$pv <- pv
  object$pe <- pe

  ## save parameters for extrapolation
  object$parex <- list(k=k,s=s,sp=sp)
  names(object$parex$k) <- pvnames

  ## chisq p-value
  if(!is.null(f <- object$fi$sphe.3)) {
    y <- sbpv1(f,get(f$psi),k=0,s=s,sp=sp)
    object$chisq <- y
  }

  ## find the best model
  aic <- sapply(object$fi,"[[","aic")
  aic0 <- min(aic)
  i <- which(aic==aic0)[1]
  model <- models[[i]]
  pv <- pv[i,]; pe <- pe[i,]
  names(pv) <- names(pe) <- pvnames
  object$best <- list(model=model,aic=aic0,pv=pv,pe=pe)

  object
}



## print
print.summary.scaleboot <- function(x,...) {

  ### raw
  a <- catpval(x$raw$pv,x$raw$pe)
  cat("\nRaw Bootstrap Probability: ",a$value,"\n")

  ## in case no fitting
  if(is.null(x$fi)) {
    cat("\nNo Model Fitting\n")
    return(invisible(x))
  }

  ## chisq
  if(!is.null(x$chisq)) {
    a <- catpval(x$chisq$pv,x$chisq$pe)
    cat("\nChisquare P-value: ",a$value)
    f <- x$fi$sphe.3
    p <- parsphere(f$par*f$mag)
    cat(" ; v=",format(p[1],digits=3),
        ", a=",format(p[2],digits=3),
        ", nu=",format(p[3],digits=3),"\n",sep="")
  }
  
  ## corrected p-values
  pval <- matrix("",nrow(x$pv),ncol(x$pv))
  dimnames(pval) <- dimnames(x$pv)
  for(i in seq(length=ncol(x$pv))) {
    a <- catpval(x$pv[,i],x$pe[,i])
    pval[,i] <- a$value
  }
  cat("\nCorrected P-values (",a$name,"):\n",sep="")
  aic <- myformat(c(pi,sbaic(x)),digits=2)[-1]
  catmat(cbind(pval,aic))

  ## best model
  cat("\nBest Model: ",x$best$model,"\n")

  invisible(x)
}


##
## scalebootv
##

summary.scalebootv <- function(object,models=attr(object,"models"),k=1:3,...) {
  for(i in seq(along=object)) object[[i]] <- summary(object[[i]],models,k=k,...)
  class(object) <- c("summary.scalebootv",class(object))
  attr(object,"models") <- models
  attr(object,"pvalues") <- paste("k",k,sep=".")
  object
}

print.summary.scalebootv <- function(x,...) {
  ## extract information
  pvalues <- attr(x,"pvalues")
  bests <- lapply(x,"[[","best")
  raws <- lapply(x,"[[","raw")
  chisqs <- lapply(x,"[[","chisq")
  nochisq <- all(sapply(chisqs,is.null))

  ## best model
  out <- matrix("",length(x),3+length(pvalues)+!nochisq,
              dimnames=list(names(x),c("raw",pvalues,if(nochisq) NULL else "chisq",
                "model","aic")))
  out[,"model"] <- format(sapply(bests,"[[","model"))
  out[,"aic"] <- format(round(sapply(bests,"[[","aic"),digits=2))
  pv <- sapply(raws,"[[","pv")
  pe <- sapply(raws,"[[","pe")
  out[,"raw"] <- catpval(pv,pe)$value
  if(!nochisq) {
    chisqs <- lapply(chisqs,function(a) if(is.null(a)) list(pv=NA,pe=NA) else a)
    pv <- sapply(chisqs,"[[","pv")
    pe <- sapply(chisqs,"[[","pe")
    out[,"chisq"] <- catpval(pv,pe)$value
  }
  for(p in pvalues) {
    pv <- sapply(bests,function(b) b$pv[[p]])
    pe <- sapply(bests,function(b) b$pe[[p]])
    out[,p] <- catpval(pv,pe)$value
  }

  cat("\nCorrected P-values (",catpval(0)$name,"):\n",sep="")
  catmat(out)

  invisible(x)
}

#######
##
## extract p-values

## general
sbpval <- function(x,...) UseMethod("sbpval")

## scaleboot
sbpval.summary.scaleboot <- function(x,sd=FALSE,best=FALSE,...) {
  if(best) {
    pv <- x$best$pv
    pe <- x$best$pe
  } else {
    pv <- x$pv
    pe <- x$pe
  }
  if(sd) {
    pval <- list(estimate=pv,sd=pe)
  } else {
    pval <- pv
  }
  pval
}

## scalebootv
sbpval.summary.scalebootv <- function(x,...) {
  lapply(x,sbpval,...)
}
