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
###
### OPTIONS
###

## model names
sbmodelnames <- function(m=1:3,one.sided=TRUE,two.sided=FALSE,rev.sided=FALSE,
                         poly,sing,poa,pob,sia,sib,sphe){
  if(missing(poly)) poly <- if(one.sided) m else 0
  if(missing(sing)) sing <- if(one.sided) m else 0
  if(missing(poa)) poa <- if(two.sided) m else 0
  if(missing(sia)) sia <- if(two.sided) m else 0
  if(missing(pob)) pob <- if(rev.sided) m else 0
  if(missing(sib)) sib <- if(rev.sided) m else 0
#  if(missing(sphe)) sphe <- if(one.sided) m else 0
  if(missing(sphe)) sphe <- 0

  make1 <- function(na,mm,min=1,max=Inf) {
    k <- mm[mm>=min & mm<=max]
    if(length(k)>0) paste(na,k,sep=".") else NULL
  }
  
  ## one sided models
  polyk <- make1("poly",poly,1)
  singk <- make1("sing",sing,3)

  ## two sided models (polynomial-difference)
  poa1k <- make1("poa1",poa,2)
  poa2k <- make1("poa2",poa,4)
  pob1k <- make1("pob1",pob,2)
  pob2k <- make1("pob2",pob,4)

  ## two sided models (singular-difference)
  sia1k <- make1("sia1",sia,4)
  sia2k <- make1("sia2",sia,5)
  sib1k <- make1("sib1",sib,4)
  sib2k <- make1("sib2",sib,5)

  ## specialized models
  sphek <- make1("sphe",sphe,3,3)

  ## output model names
  c(polyk,singk,poa1k,poa2k,pob1k,pob2k,sia1k,sia2k,sib1k,sib2k,sphek)
}

## default options (local variable)

.onLoad <- function(libname, pkgname) {
  op <- list(
  ## sbfit
  models = sbmodelnames(), # default models for fitting
  control = list(reltol=1e-14),  # for optim
  method = NULL,  # for optim
  ## poly model
  mag.poly = c(1,0.1,0.01,0.001), # mag factor for par in poly model
  omg.poly = c(0.1,1,1,1), # weights
  ## sing model
  mag1.sing = 0.1, # mag factor for singularity parameter
  omg1.sing = 1.0, # weight 
  ## sphe model
  mag.sphe = c(1.0,1.0,0.01), # mag factor for spherical model
  omg.sphe = c(0.1,1.0,1.0), # weights
  chisq.sphe = FALSE, # use chisq for pvalues instead of extrapolation for k >= 2
  ## dif parms of poa and pob models
  mag1.poa = c(1,0.1),  # mag factor for dif parm
  low1.poa = c(1,-3), # lower limits of dif parm
  upp1.poa = c(7,3), # upper limits of dif parm
  trg1.poa = c(10,0), # targets of dif parm
  omg1.poa = c(1,1), # weights of dif parm
  ## print parameters
  percent = TRUE, # print p-values in percent
  digits.pval = 2, # significant digits for pvalue
  digits.coef = 4, # significant digits for coefficients
  ## averaging models
  th.aicw = 1e-4,  # threshold for akaike weights

  ## sbconf; only experimental...
  probs0=c(0.001,0.01,0.1,0.9,0.99,0.999), # initial grid
  prob0=0.5, # for starting value
  tol.mono=0.1,  # for monotonicity checking
  tol.conv=0.01, # for convergence
  tol.relconv=0.5, # tolerance with respect to s.e. 
  max.loop=100, # max iterations

  ## for debug
  debug=FALSE
  )
  options(scaleboot=op)
}


## to show and change options
sboptions <- function(x,value) {
  op <- getOption("scaleboot")
  if(missing(x)) return(op)
  if(length(x)==1) y <- op[[x]]
  else y <- op[x]
  if(!missing(value)) {
    op[[x]] <- value
    options(scaleboot=op)
  }
  y
}

######################################################################
### EOF
