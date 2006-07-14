## makedata.R

date()
library(snow)
library(pvclust)
library(scaleboot)
data(lung)
cl <- makeCluster(40,"MPI")
sa <- 10^seq(-2,2,length=13) # wider range of scales than pvclust default
lung73.pvclust <- parPvclust(cl,lung,r=1/sa,nboot=10000,weight=T) 
lung73.sb <- sbfit(lung73.pvclust,cluster=cl) # model fitting
save(lung73.pvclust,lung73.sb,file="lung73.rda")


date()
mam15.mt <- read.mt("mam15.mt")
mam15.ass <- read.ass("mam15.ass")
set.seed(100)
mam15.relltest <- relltest(mam15.mt,nb=100000,ass=mam15.ass)
save(mam15.mt,mam15.ass,mam15.relltest,file="mam15.rda")

date()
quit(save="no")
