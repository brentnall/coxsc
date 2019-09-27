## Summary of results from simulations

##load functions
source("pifuncs.R");  dyn.load("../libs/coxsc.so")

#to format numbers to text for table
fn.format<-function(ind, ndigit=3)
  {
            format(round(ind,ndigit), nsmall=ndigit)
          }


myn<-1000

my2par<-c(0.2, 0.7, 0.05, 0.05) #tau, theta, cen before, cen after

insp.1<-list(myn=myn, mypi=0.25, mypar=list(tau=my2par[1], theta=my2par[2], base=fn.base(my2par)), mycen=fn.cen(my2par))

my2par<-c(0.2, 0.7, 0.1, 0.9) #tau, theta, cen before, cen after

insp.2<-list(myn=myn, mypi=0.25, mypar=list(tau=my2par[1], theta=my2par[2], base=fn.base(my2par)), mycen=fn.cen(my2par))

my2par<-c(2, 0.7, 0.05, 0.05) #tau, theta, cen before, cen after

insp.1b<-list(myn=myn, mypi=0.25, mypar=list(tau=my2par[1], theta=my2par[2], base=fn.base(c(1,my2par[2:4]))), mycen=fn.cen(c(1,my2par[2:4])))

my2par<-c(2, 0.7, 0.1, 0.9) #tau, theta, cen before, cen after

insp.2b<-list(myn=myn, mypi=0.25, mypar=list(tau=my2par[1], theta=my2par[2], base=fn.base(c(1,my2par[2:4]))), mycen=fn.cen(c(1,my2par[2:4])))
         
my2par<-c(0.2, 1.4, 0.05, 0.05) #tau, theta, cen before, cen after

insp.1c<-list(myn=myn, mypi=0.25, mypar=list(tau=my2par[1], theta=my2par[2], base=fn.base(c(1,my2par[c(1,3,4)]))), mycen=fn.cen(c(1,my2par[c(1,3,4)])))

my2par<-c(0.2, 1.4, 0.1, 0.9) #tau, theta, cen before, cen after

insp.2c<-list(myn=myn, mypi=0.25, mypar=list(tau=my2par[1], theta=my2par[2], base=fn.base(c(1,my2par[c(1,3,4)]))), mycen=fn.cen(c(1,my2par[c(1,3,4)])))
    
my2par<-c(5, 1.4, 0.05, 0.05) #tau, theta, cen before, cen after

insp.1d<-list(myn=myn, mypi=0.25, mypar=list(tau=my2par[1], theta=my2par[2], base=fn.base(c(1,1,my2par[c(3,4)]))), mycen=fn.cen(c(1,1,my2par[c(3,4)])))

my2par<-c(5, 1.4, 0.1, 0.9) #tau, theta, cen before, cen after

insp.2d<-list(myn=myn, mypi=0.25, mypar=list(tau=my2par[1], theta=my2par[2], base=fn.base(c(1,1,my2par[c(3,4)]))), mycen=fn.cen(c(1,1,my2par[c(3,4)])))

insp<-list(insp.1, insp.2, insp.1b, insp.2b, insp.1c, insp.2c, insp.1d, insp.2d)


###############
##results
###############
load("simres_scen22.Rdata"); load("simresFL_scen2.Rdata"); load("simresFL_scen22.Rdata"); load("simresFL_scen222.Rdata")

nreps<-length(myout1[[1]])

myout2<-lapply(myout1,function(ind)matrix(unlist(ind), ncol=nreps))

load("simres_scen22.Rdata")

myout3<-lapply(myout1,function(ind)matrix(unlist(ind), ncol=nreps))

#some plots
par(mfrow=c(1,2))

boxplot(lapply(myout2, function(ind) exp(ind[2,])))

boxplot(lapply(myout3, function(ind) exp(ind[2,])))

boxplot(lapply(myout2, function(ind) (ind[3,])))

boxplot(lapply(myout3, function(ind) (ind[3,])))

##percent bias

mymnbias<-sapply(1:length(insp), function(idx) rowMeans ((myout2[[idx]][1:2,] - log(unlist(insp[[idx]]$mypar[1:2])))/log(unlist(insp[[idx]]$mypar[1:2]))))*100

mymnbiaspl<-sapply(1:length(insp), function(idx) rowMeans ((myout3[[idx]][1:2,] - log(unlist(insp[[idx]]$mypar[1:2])))/log(unlist(insp[[idx]]$mypar[1:2]))))*100

#log variance

mylogvar<-sapply(myout2, function(ind)apply(ind[1:2,],1,var))

mylogvarpl<-sapply(myout3, function(ind)apply(ind[1:2,],1,var))

##u
mymnu<-sapply(1:length(insp), function(idx) mean(myout2[[idx]][3,]))*100

mymnupl<-sapply(1:length(insp), function(idx) mean(myout3[[idx]][3,]))*100

mysdu<-sapply(1:length(insp), function(idx) sd(myout2[[idx]][3,]*100))

mysdupl<-sapply(1:length(insp), function(idx) sd(myout3[[idx]][3,]*100))


##table
mytab.bias<-rbind(mymnbias[2,], mymnbiaspl[2,], mymnbias[1,], mymnbiaspl[1,])

mytab.var<-rbind(mylogvar[2,], mylogvarpl[2,], mylogvar[1,], mylogvarpl[1,])*100

mytab.u<-rbind(mymnu, mymnupl, mysdu, mysdupl)

mytab<-rbind(mytab.bias, mytab.var, mytab.u)

fn.format<-function(ind, ndigit=3)
  {
        format(round(ind,ndigit), nsmall=ndigit)
      }


mytab<-fn.format(mytab)

library("xtable")

xtable(mytab)

##########bootstrap results
load("simresBOOT_scen22.Rdata")

nreps<-length(myout1[[1]])

myout2<-lapply(myout1,function(ind)matrix(unlist(ind), ncol=nreps))

## fits, real number, bootstrap variance and percentile ci's

obsvar<-sapply(myout2, function(ind) apply(ind[1:2,],1,var)) #observed variance

bsvar<-sapply(myout2, function(ind) apply(ind[6:7,],1,median)) #median bootstrap variance

bsvar<-sapply(myout2, function(ind) apply(ind[6:7,],1,mean)) #mean bootstrap variance

bsvarsd<-sapply(myout2, function(ind) apply(ind[6:7,],1,sd)) #sd bootstrap 

mytab<-rbind(rbind(obsvar[2,], bsvar[2,], bsvarsd[2,])*100,rbind(obsvar[1,], bsvar[1,], bsvarsd[1,])*100)

mytab<-fn.format(mytab)

xtable(mytab)

