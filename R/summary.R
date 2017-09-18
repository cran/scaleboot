##
##  scaleboot: R package for multiscale bootstrap
##  Copyright (C) 2006-2008 Hidetoshi Shimodaira
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
            pvalues = attr(x,"pvalues"),
            lambda = attr(x,"lambda")
            )
##
## scaleboot
##

summary.scaleboot <- function(object,models=names(object$fi),
                              k=3,s=1,sp=-1,
                              type=c("Frequentist","Bayesian"),...) {
## note: ... is not passed to any further, but only to avoid error
##                 "S3 generic/method consistency ... WARNING"

  ## option
  op <- sboptions()

  ## object type
  class(object) <- c("summary.scaleboot",class(object))

  ## p-values
  pvnames <- paste("k",k,sep=".")
  if(!is.numeric(type)) {
    type <- match.arg(type)
    lambda <- switch(type, Bayesian=0, Frequentist=1)
  } else lambda <- type

  ## save parameters for extrapolation
  object$parex <- list(k=k,s=s,sp=sp,lambda=lambda)
  names(object$parex$k) <- pvnames

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
    object$best <- list(model="raw",aic=0,pv=pv,pe=pe,dv=0,de=0)
    object$average <- list(model="raw",w=structure(1,names="raw"),pv=pv,pe=pe,dv=0,de=0)
    return(object)
  }

  ## corrected p-values and beta0
  pv <- pe <- matrix(NA,length(models),length(k),
                     dimnames=list(models,pvnames))
  dv <- de <- structure(rep(0,length(models)),names=models) # beta0 (value and sd)
  for(i in seq(along=models)) {
    m <- models[[i]]
    f <- object$fi[[m]]
    psi <- sbpsiget(f$model)
    for(j in seq(along=k)) {
      y <- sbpv1(f,psi,k=k[[j]],s=s,sp=sp,lambda=lambda)
      pv[i,j] <- y$pv
      pe[i,j] <- y$pe
    }
    dv[i] <- f$par[1]*f$mag[1]
    de[i] <- sqrt(f$var[1,1])*f$mag[1]
  }
  object$pv <- pv # p-value
  object$pe <- pe # (sd)
  object$dv <- dv # beta0
  object$de <- de # (sd)

  ## chisq p-value
  if(!is.null(f <- object$fi$sphe.3)) {
    y <- sbpv1(f,sbpsiget(f$model),k=0,s=s,sp=sp)
    object$chisq <- y
  }

  ## find the best model
  aic <- sapply(object$fi,"[[","aic")
  aic0 <- min(aic)
  i <- which(aic==aic0)[1]
  model <- models[[i]]
  pvbest <- pv[i,]; pebest <- pe[i,]
  names(pvbest) <- names(pebest) <- pvnames
  object$best <- list(model=model,aic=aic0,pv=pvbest,pe=pebest,dv=dv[[i]],de=de[[i]])

  ## average by akaike weights
  w <- exp(-(aic-aic0)/2) # akaike weights
  w <- w/sum(w)
  u <- w>op$th.aicw # ignore small values
  w <- w[u]/sum(w[u])
  pvave <- apply(w*pv[u,,drop=F],2,sum)
  peave <- apply(w*pe[u,,drop=F],2,sum)
  dvave <- sum(w*dv[u])
  deave <- sum(w*de[u])
  object$average <- list(w=w,pv=pvave,pe=peave,dv=dvave,de=deave)

  object
}

### corrected p-value 1
## fit : output of optims (includes var, mag)
## psi : function(beta,s,k)
## k : degree for corrected p-value (default: k=1)
## s : sigma^2 for corrected p-value (default: s=1)
## sp : sigma^2 for prediction (default: sp=-1)
## lambda : mixing bayes (lambda=0) and freq (lambda=1)

sbpv1 <- function(fit,psi,k=1,s=1,sp=-1,lambda=0) {
  pval <- function(par) pnorm(-psi(fit$mag*par,s,k=k,sp=sp,lambda=lambda))
  pv <- pval(fit$par)
  h <- nderiv(pval,fit$par)
  pe <- sqrtx(h %*% fit$var %*% h)
  list(pv=pv,pe=pe)
}

## print
print.summary.scaleboot <- function(x,sort.by=c("aic","none"),verbose=FALSE,...) {
  ## verbose
  if(verbose) print.scaleboot(x,sort.by=sort.by,...)
  
  ### raw
  a <- catpval(x$raw$pv,x$raw$pe)
  if(is.null(x$raw$s)) x$raw$s <- NA
  cat("\nRaw Bootstrap Probability (scale=",round(x$raw$s,3),
      ") : ",a$value,"\n",sep="")

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

  ## corrected p-values (for models, and also for the best and average)
  pvs <- rbind(x$best$pv,x$average$pv,x$pv)
  pes <- rbind(x$best$pe,x$average$pe,x$pe)
  rownames(pvs)[1:2] <- rownames(pes)[1:2] <- c("best","average")

  ## beta0
  dvs <- c(x$best$dv,x$average$dv,x$dv)
  des <- c(x$best$de,x$average$de,x$de)
  names(dvs)[1:2] <- names(des)[1:2] <- c("best","average")

  ## prepare table
  pval <- matrix("",nrow(pvs),ncol(pvs))
  dimnames(pval) <- dimnames(pvs)
  
  for(i in seq(length=ncol(pvs))) {
    a <- catpval(pvs[,i],pes[,i],lambda=x$parex$lambda)
    pval[,i] <- a$value
  }
  cat("\nCorrected P-values for Models (",a$name,",",a$lambda,"):\n",sep="")
  pvalbest <- pval[1:2,,drop=F] # for best and average
  pval <- pval[-(1:2),,drop=F] # for models

  ## beta0 table
  beta0 <- myformat(c(pi,dvs),c(pi,des),digits=2)[-1]
  beta0best <- beta0[1:2]
  beta0 <- beta0[-(1:2)]

  ## aic and akaike weights
  aicval <- sbaic(x)
  aic <- myformat(c(pi,aicval),digits=2)[-1]
  weight <- aic; weight[] <- ""
  a <- catpval(x$average$w)$value
  weight[names(a)] <- a

  ## sort
  tab <- cbind(pval,beta0,aic,weight) # to be catmat
  sort.by <- match.arg(sort.by)
  j <- switch(sort.by,
              none=1:length(aicval),
              aic=order(aicval))
  tabj <- tab[j,]

  ## print the table
  catmat(tabj)

  ## best model
  cat("\nBest Model: ",x$best$model,"\n")

  ## the bottom line
  cat("\nCorrected P-values by the Best Model and by Akaike Weights Averaging:\n")
  beta0 <- beta0best
  catmat(cbind(pvalbest,beta0))
  cat("\n")

  invisible(x)
}


##
## scalebootv
##

summary.scalebootv <- function(object,models=attr(object,"models"),k=3,
                               type="Frequentist",...) {
  for(i in seq(along=object)) object[[i]] <- summary(object[[i]],models,k=k,type=type,...)
  class(object) <- c("summary.scalebootv",class(object))
  attr(object,"models") <- models
  attr(object,"pvalues") <- paste("k",k,sep=".")
  attr(object,"lambda") <- object[[1]]$parex$lambda
  object
}


selectpv <- function(x,select) {
  models <- attr(x,"models")
  select <- match.arg(select,c("best","average",models))
  if(select=="best") {
    pvpe <- lapply(x,"[[","best")
    model <- format(sapply(pvpe,"[[","model"))
    aic <- format(round(sapply(pvpe,"[[","aic"),digits=2))
    outaic <- cbind(model,aic)
    selna <- "the Best Model"
  } else if(select=="average") {
    pvpe <- lapply(x,function(s)
       list(model=s$best$model,weight=s$average$w[s$best$model],
            pv=s$average$pv,pe=s$average$pe,
            dv=s$average$dv,de=s$average$de))
    model <- format(sapply(pvpe,"[[","model"))
    weight <- catpval(sapply(pvpe,"[[","weight"))$value
    outaic <- cbind(model,weight)
    selna <- "Akaike Weights Averaging"
  } else {
    pvpe <- lapply(x, function(s)
      if(!is.null(s$fi)) list(model=select,aic=s$fi[[select]]$aic,
                              pv=s$pv[select,],pe=s$pe[select,],
                              dv=s$dv[select],de=s$de[select])
      else s$best)
    model <- format(sapply(pvpe,"[[","model"))
    aic <- format(round(sapply(pvpe,"[[","aic"),digits=2))
    outaic <- cbind(model,aic)
    selna <- select
  }
  
  list(pvpe=pvpe,select=select,name=selna,outaic=outaic)
}


print.summary.scalebootv <- function(x,select="average",sort.by=NULL,nochisq=TRUE,...) {
  ## extract information
  pvalues <- attr(x,"pvalues")
  lambda <- attr(x,"lambda")
  raws <- lapply(x,"[[","raw")
  if(!nochisq) {
    chisqs <- lapply(x,"[[","chisq")
    nochisq <- all(sapply(chisqs,is.null))
  }

  ## prepare table containers for p-values
  out <- matrix("",length(x),2+length(pvalues)+!nochisq,
              dimnames=list(names(x),c("raw",pvalues,if(nochisq) NULL else "chisq","beta0")))

  outval <- matrix(0,length(x),1+length(pvalues)+!nochisq,
                   dimnames=list(names(x),c("raw",pvalues,if(nochisq) NULL else "chisq"))
                   ) # numerical values for sorting

  ## which p-values to be printed?
  selpv <- selectpv(x,select)
  out <- cbind(out,selpv$outaic)

  ## fill-in values
  pv <- sapply(raws,"[[","pv")
  pe <- sapply(raws,"[[","pe")
  out[,"raw"] <- catpval(pv,pe)$value
  outval[,"raw"] <- pv
  if(!nochisq) {
    chisqs <- lapply(chisqs,function(a) if(is.null(a)) list(pv=NA,pe=NA) else a)
    pv <- sapply(chisqs,"[[","pv")
    pe <- sapply(chisqs,"[[","pe")
    out[,"chisq"] <- catpval(pv,pe)$value
    outval[,"chisq"] <- pv
  }
  for(p in pvalues) {
    pv <- sapply(selpv$pvpe,function(b) b$pv[[p]])
    pe <- sapply(selpv$pvpe,function(b) b$pe[[p]])
    a <- catpval(pv,pe,lambda=lambda)
    out[,p] <- a$value
    outval[,p] <- pv
  }
  dv <- sapply(selpv$pvpe,"[[","dv")
  de <- sapply(selpv$pvpe,"[[","de")
  de[de>10] <- NA
  out[,"beta0"] <- myformat(c(pi,dv),c(pi,de),digits=2)[-1]

  ## sort and print
  cat("\nCorrected P-values by ",selpv$name," (",a$name,",",a$lambda,"):\n",sep="")
  if(!is.null(sort.by) && sort.by!="none") {
    j <- order(-outval[,sort.by])
    catmat(out[j,])
  } else catmat(out)

  invisible(x)
}

#######
##
## extract p-values

## general
sbpval <- function(x,...) UseMethod("sbpval")

## scaleboot
sbpval.summary.scaleboot <- function(x,sd=FALSE,
                                     select=c("average","best","all"),...) {
  select <- match.arg(select)
  y <- switch(select,
              average=x$average,
              best=x$best,
              all=x)
  pv <- y$pv
  pe <- y$pe

  if(sd) {
    pval <- list(estimate=pv,sd=pe)
  } else {
    pval <- pv
  }
  pval
}

## scalebootv
sbpval.summary.scalebootv <- function(x,...) {
  sapply(x,sbpval,...)
}
