## Fit a model using partial likelihood

source("pifuncs.R"); dyn.load("../libs/coxsc.so")

insp1<-fn.simscen(1000)

myout<- fn.plsim(1, insp1[[1]])

