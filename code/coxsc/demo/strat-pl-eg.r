## Example: fit stratified model

## Load functions
source("fitfns.r");    dyn.load("../libs/coxsc.so"); library("boot"); library("survival")

##number replicates
nreps<-1

##generate data stratified
myfiles<- sapply(1:nreps, function(idx) paste("bigsims_", idx))
temp<-sapply(1:nreps, function(idx) fn.genbigsim(idx, infile=myfiles[idx]))

##no ci
myout1<-sapply(1:nreps, function(idy) fn.main.stratified(myfiles[idy], doci=FALSE))
save(file="bigsimres.Rdata",myout1)

## ci
myout2<-sapply(1:nreps, function(idy) fn.main.stratified(myfiles[idy], doci=TRUE))
save(file="bigsimres2.Rdata",myout2)

