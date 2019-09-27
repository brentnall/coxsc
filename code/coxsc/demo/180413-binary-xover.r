## Two arm trial with selective crossover, binary models
## Analysis reported in the paper
## Adam Brentnall
## 13th April 2018

########################
## 1. Functions
########################

fn.solvepi<-function(p.n, p.ne, p.nc, p.w, p.us)
  {
    mya<- (1 - p.w) * (p.nc / p.n - 1) / p.n
    myb<-1 + p.us * (1 - p.w) / p.n - (p.ne * p.w + p.nc) / p.n
    myc<--p.us 

    ##polynomial solve
    u0<- c( (-myb + sqrt(myb^2 - 4*mya*myc) ) / (2*mya), (-myb - sqrt(myb^2 - 4*mya*myc) ) / (2*mya))

    u0[1]/p.n
    
  }

fn.fitu1<-function(p.n, p.ne, p.nc, p.w, p.u0)
  {

     p.u1<-p.u0 - (p.ne * (p.u0 / p.n) * p.w) / (1 - (p.u0 / p.n) * (1 - p.w)) - p.nc * p.u0 / p.n
     p.u1

  }


## Likelihood
### calc pi as part of lik routine each time
fn.lik2<-function(mypar,mydta)
  {
  
     p.omega<-exp(mypar[3]); p.theta<-exp(mypar[4]); p.base <- exp(mypar[1]); p.base2 <- exp(mypar[2])

     ##myparmax<-10 # limit range of exploration
     ##p.omega<-myparmax / (1+p.omega); p.theta<-myparmax / (1+p.omega)
     ##p.base <- 1 / (1+p.base); p.base2 <- 1 / (1+p.base2)
     
     t.n0<-mydta[2:3]; t.n1<-mydta[4:6]; t.y0<-mydta[7:8]; t.y1<-mydta[9:11]

     p.pi<-fn.solvepi(t.n0[2], t.y0[2], t.n0[2] - t.y0[2] - t.n1[2] -t.n1[3], p.omega, t.n1[2])

      
##     print(p.pi)
     
     ## number switchers in letrozole at switch, updates    
    p.pi1<-fn.fitu1(t.n0[1], t.y0[1], t.n0[1] - t.n1[1] - t.y0[1], p.omega, p.pi * t.n0[1]) / t.n1[1]

    t.pL0<-p.pi* p.base*p.theta*p.omega + (1-p.pi) * p.base * p.theta
    t.pL1<-p.pi1* p.base2*p.theta*p.omega + (1-p.pi1) * p.base2 * p.theta
    
    t.pT0<-p.pi* p.base*p.omega + (1-p.pi) * p.base 
    t.pT1<-c(p.base2*p.omega*p.theta, p.base2) #switch don't       
    
    mylikL0 <- t.y0[1] * log(t.pL0) +  (t.n0[1] - t.y0[1]) * log(1 - t.pL0)
    mylikL1 <- t.y1[1] * log(t.pL1) +  (t.n1[1] - t.y1[1]) * log(1 - t.pL1)
    
    mylikT0 <- t.y0[2] * log(t.pT0) +  (t.n0[2] - t.y0[2]) * log(1 - t.pT0)
    mylikT1x <- t.y1[2] * log(t.pT1[1]) +  (t.n1[2] - t.y1[2]) * log(1 - t.pT1[1])
    mylikT1r <- t.y1[3] * log(t.pT1[2]) +  (t.n1[3] - t.y1[3]) * log(1 - t.pT1[2])
    
    mylik<-c(mylikL0, mylikL1, mylikT0, mylikT1x, mylikT1r)

    -sum(mylik)

  }

## likelihood used for profile - theta in the data part
fn.likprof2<-function(mypar, mydta)
  {
    mypar1<-c(mypar, mydta[1]); mydta1<-mydta[2:length(mydta)]
    fn.lik2(mypar1, mydta1)

 }


## fit model with pi in lik calc
fn.fitme2<-function(initpar, initdta)
  {

     p.pi0<-initdta[1]; t.n0<-initdta[2:3]; t.n1<-initdta[4:6]; t.y0<-initdta[7:8]; t.y1<-initdta[9:11]

     mydta<-c(NA, t.n0, t.n1, t.y0, t.y1); mypar<-initpar

     t.opt<-optim(par=mypar, fn=fn.lik2, mydta=mydta, method="BFGS")

##     print(t.opt)
     
     c(t.opt$value, exp(t.opt$par))
    
  }



## fit model, fix theta though. used for profile lik ci
fn.fitprofile2<-function(initpar, initdta, niter=10, mytol = 0.000001)
  {

     p.pi0<-initdta[1]; t.n0<-initdta[2:3]; t.n1<-initdta[4:6]; t.y0<-initdta[7:8]; t.y1<-initdta[9:11]

     p.pi<-fn.solvepi(t.n0[2], t.y0[2], t.n0[2] - t.y0[2] - t.n1[2] -t.n1[3], exp(initpar[3]), t.n1[2])

     mydta<-c(initpar[4], p.pi, t.n0, t.n1, t.y0, t.y1)

     t.opt<-optim(par=initpar[1:3], fn=fn.likprof2, mydta=mydta, method="BFGS")

     c(t.opt$value, exp(t.opt$par))
    
  }

##used by uniroot to find confindence intervals
fn.profroot2<-function(inpar, myML, mychi, mypar, mydta)
  {
    mypfit<-fn.fitprofile2(c(mypar[1:3], mypar[4]+inpar), mydta)
  - myML + mypfit[1] - mychi
  }


## Profile likelihood CI
fn.profCI2<-function(myfit, mydta, myalpha=0.05)
{
  mychi<-qchisq(1-myalpha, 1)/2

  myCIU<-uniroot(fn.profroot2, interval=c(0,0.5), myML=myfit[1], mychi=mychi, mypar=log(myfit[2:5]), mydta=mydta)
  myCIL<-uniroot(fn.profroot2, interval=c(-0.5,0), myML=myfit[1], mychi=mychi, mypar=log(myfit[2:5]), mydta=mydta)
  myCI<-exp(c(myCIL$root, myCIU$root) + log(myfit[5]))

  myrng<-seq(myCIL$root-0.1, myCIU$root+0.1, by = (myCIU$root- myCIL$root + 0.2)/20)
  myPL<-sapply(myrng, function(idx) fn.fitprofile2(c(log(myfit[2:4]), log(myfit[5])+idx), mydta))

  plot(exp(log(myfit[5])+myrng), -myPL[1,], xlab="Treatment effect", ylab="Profile log likelihood")
  points(c(myCI[1], myCI[2], myfit[5]), c(-myfit[1]-mychi, -myfit[1]-mychi, -myfit[1]), col=2)
  lines(c(myCI[1], myCI[1]), c(-10000,-myfit[1]-mychi), lty=2)
  lines(c(myCI[2], myCI[2]), c(-10000,-myfit[1]-mychi),lty=2)
  lines(c(myfit[5], myfit[5]), c(-10000,-myfit[1]), lwd=2)

  myCI

}


############ TWO TREATMENT EFFECTS
## Likelihood, time-dept treatment
### calc pi as part of lik routine each time
fn.lik3<-function(mypar,mydta)
  {
  
     p.omega<-exp(mypar[3]); p.theta<-exp(mypar[4]); p.base <- exp(mypar[1]); p.base2 <- exp(mypar[2]); p.theta2<-exp(mypar[5]);

     ##myparmax<-10 # limit range of exploration
     ##p.omega<-myparmax / (1+p.omega); p.theta<-myparmax / (1+p.omega)
     ##p.base <- 1 / (1+p.base); p.base2 <- 1 / (1+p.base2)
     
     t.n0<-mydta[2:3]; t.n1<-mydta[4:6]; t.y0<-mydta[7:8]; t.y1<-mydta[9:11]

     p.pi<-fn.solvepi(t.n0[2], t.y0[2], t.n0[2] - t.y0[2] - t.n1[2] -t.n1[3], p.omega, t.n1[2])

##     print(p.pi)
     
     ## number switchers in letrozole at switch, updates    
    p.pi1<-fn.fitu1(t.n0[1], t.y0[1], t.n0[1] - t.n1[1] - t.y0[1], p.omega, p.pi * t.n0[1]) / t.n1[1]

    t.pL0<-p.pi* p.base*p.theta*p.omega + (1-p.pi) * p.base * p.theta
    t.pL1<-p.pi1* p.base2*p.theta2*p.omega + (1-p.pi1) * p.base2 * p.theta2
    
    t.pT0<-p.pi* p.base*p.omega + (1-p.pi) * p.base 
    t.pT1<-c(p.base2*p.omega*p.theta2, p.base2) #switch don't       
    
    mylikL0 <- t.y0[1] * log(t.pL0) +  (t.n0[1] - t.y0[1]) * log(1 - t.pL0)
    mylikL1 <- t.y1[1] * log(t.pL1) +  (t.n1[1] - t.y1[1]) * log(1 - t.pL1)
    
    mylikT0 <- t.y0[2] * log(t.pT0) +  (t.n0[2] - t.y0[2]) * log(1 - t.pT0)
    mylikT1x <- t.y1[2] * log(t.pT1[1]) +  (t.n1[2] - t.y1[2]) * log(1 - t.pT1[1])
    mylikT1r <- t.y1[3] * log(t.pT1[2]) +  (t.n1[3] - t.y1[3]) * log(1 - t.pT1[2])
    
    mylik<-c(mylikL0, mylikL1, mylikT0, mylikT1x, mylikT1r)

    -sum(mylik)

  }


## fit model with pi in lik calc
fn.fitme3<-function(initpar, initdta)
  {

     p.pi0<-initdta[1]; t.n0<-initdta[2:3]; t.n1<-initdta[4:6]; t.y0<-initdta[7:8]; t.y1<-initdta[9:11]

     mydta<-c(NA, t.n0, t.n1, t.y0, t.y1); mypar<-initpar

     t.opt<-optim(par=mypar, fn=fn.lik3, mydta=mydta, method="BFGS")

##     print(t.opt)
     
     c(t.opt$value, exp(t.opt$par))
    
  }

## likelihood used for profile - theta in the data part
fn.likprof3<-function(mypar, mydta, whichfix=0)
  {
      if(whichfix==0){ ##fix period 0
          mypar1<-c(mypar[1:3], mydta[1], mypar[4]); mydta1<-mydta[2:length(mydta)]
      } else { #fix period 0 
          mypar1<-c(mypar, mydta[1]); mydta1<-mydta[2:length(mydta)]
      }
          
    fn.lik3(mypar1, mydta1)

 }



## fit model, fix theta0 though. used for profile lik ci
fn.fitprofile3<-function(initpar, initdta, niter=10, mytol = 0.000001, whichfix=0)
  {

     p.pi0<-initdta[1]; t.n0<-initdta[2:3]; t.n1<-initdta[4:6]; t.y0<-initdta[7:8]; t.y1<-initdta[9:11]

     p.pi<-fn.solvepi(t.n0[2], t.y0[2], t.n0[2] - t.y0[2] - t.n1[2] -t.n1[3], exp(initpar[3]), t.n1[2])


     if(whichfix==0){
         mydta<-c(initpar[4], p.pi, t.n0, t.n1, t.y0, t.y1)
         t.opt<-optim(par=c(initpar[1:3], initpar[5]), fn=fn.likprof3, mydta=mydta, method="BFGS", whichfix=whichfix)
     } else {
         mydta<-c(initpar[5], p.pi, t.n0, t.n1, t.y0, t.y1)
         t.opt<-optim(par=c(initpar[1:4]), fn=fn.likprof3, mydta=mydta, method="BFGS", whichfix=whichfix)
     }

     c(t.opt$value, exp(t.opt$par))
    
  }

##used by uniroot to find confindence intervals
fn.profroot3<-function(inpar, myML, mychi, mypar, mydta, whichfix=0)
  {
      if(whichfix==0){
          mypfit<-fn.fitprofile3(c(mypar[1:3], mypar[4]+inpar, mypar[5]), mydta, whichfix=whichfix)
      } else {
          mypfit<-fn.fitprofile3(c(mypar[1:4], mypar[5]+inpar), mydta, whichfix=whichfix)
      }
  - myML + mypfit[1] - mychi
  }


## Profile likelihood CI
fn.profCI3<-function(myfit, mydta, myalpha=0.05, whichfix=0, plotrange=c(0.5,1.2), yplotrange=5)
{
#myfit3<-fn.fitme3(c(myparinit, myparinit[4]), mydta)
#-2*(myfit3[1]-myfit2[1])
##debug
##myfit<-myfit3 ; myalpha<-0.05; whichfix<-1

    mychi<-qchisq(1-myalpha, 1)/2

  myCIU<-uniroot(fn.profroot3, interval=c(0,0.5), myML=myfit[1], mychi=mychi, mypar=log(myfit[2:6]), mydta=mydta, whichfix=whichfix)
  myCIL<-uniroot(fn.profroot3, interval=c(-0.5,0), myML=myfit[1], mychi=mychi, mypar=log(myfit[2:6]), mydta=mydta, whichfix=whichfix)




  myrng<-seq(myCIL$root-0.1, myCIU$root+0.1, by = (myCIU$root- myCIL$root + 0.2)/50)

if(whichfix==0){
    myCI<-exp(c(myCIL$root, myCIU$root) + log(myfit[5]))
    myPL<-sapply(myrng, function(idx) fn.fitprofile3(c(log(myfit[2:4]), log(myfit[5])+idx, log(myfit[6])), mydta, whichfix=whichfix))
    plot(exp(log(myfit[5])+myrng), -myPL[1,], xlab="Treatment effect", ylab="Profile log likelihood", lty=1, xlim=plotrange, type="l", log="x", ylim=c(-myfit[1]-yplotrange, -myfit[1]))
    points(c(myCI[1], myCI[2], myfit[5]), c(-myfit[1]-mychi, -myfit[1]-mychi, -myfit[1]), col=2)
    lines(c(myCI[1], myCI[1]), c(-10000,-myfit[1]-mychi), lty=2)
    lines(c(myCI[2], myCI[2]), c(-10000,-myfit[1]-mychi),lty=2)
    lines(c(myfit[5], myfit[5]), c(-10000,-myfit[1]), lwd=2)
} else {
        myCI<-exp(c(myCIL$root, myCIU$root) + log(myfit[6]))
     myPL<-sapply(myrng, function(idx) fn.fitprofile3(c(log(myfit[2:5]), log(myfit[6])+idx), mydta, whichfix=1))
    plot(exp(log(myfit[6])+myrng), -myPL[1,], xlab="Treatment effect", ylab="Profile log likelihood", lty=1, xlim=plotrange, type="l", log="x",  ylim=c(-myfit[1]-yplotrange, -myfit[1]))
    points(c(myCI[1], myCI[2], myfit[6]), c(-myfit[1]-mychi, -myfit[1]-mychi, -myfit[1]), col=2)
    lines(c(myCI[1], myCI[1]), c(-10000,-myfit[1]-mychi), lty=2)
    lines(c(myCI[2], myCI[2]), c(-10000,-myfit[1]-mychi),lty=2)
    lines(c(myfit[6], myfit[6]), c(-10000,-myfit[1]), lwd=2)
}    

  myCI

}



########################################################
##model with different insistor effect
fn.lik4<-function(mypar,mydta)
  {
  
     p.omega<-exp(mypar[3]); p.theta<-exp(mypar[4]); p.base <- exp(mypar[1]); p.base2 <- exp(mypar[2]); p.omega2<-exp(mypar[5]);

     ##myparmax<-10 # limit range of exploration
     ##p.omega<-myparmax / (1+p.omega); p.theta<-myparmax / (1+p.omega)
     ##p.base <- 1 / (1+p.base); p.base2 <- 1 / (1+p.base2)
     
     t.n0<-mydta[2:3]; t.n1<-mydta[4:6]; t.y0<-mydta[7:8]; t.y1<-mydta[9:11]

     p.pi<-fn.solvepi(t.n0[2], t.y0[2], t.n0[2] - t.y0[2] - t.n1[2] -t.n1[3], p.omega, t.n1[2])

##     print(p.pi)
     
     ## number switchers in letrozole at switch, updates    
    p.pi1<-fn.fitu1(t.n0[1], t.y0[1], t.n0[1] - t.n1[1] - t.y0[1], p.omega, p.pi * t.n0[1]) / t.n1[1]

    t.pL0<-p.pi* p.base*p.theta*p.omega + (1-p.pi) * p.base * p.theta
    t.pL1<-p.pi1* p.base2*p.theta*p.omega2 + (1-p.pi1) * p.base2 * p.theta
    
    t.pT0<-p.pi* p.base*p.omega + (1-p.pi) * p.base 
    t.pT1<-c(p.base2*p.omega2*p.theta, p.base2) #switch don't       
    
    mylikL0 <- t.y0[1] * log(t.pL0) +  (t.n0[1] - t.y0[1]) * log(1 - t.pL0)
    mylikL1 <- t.y1[1] * log(t.pL1) +  (t.n1[1] - t.y1[1]) * log(1 - t.pL1)
    
    mylikT0 <- t.y0[2] * log(t.pT0) +  (t.n0[2] - t.y0[2]) * log(1 - t.pT0)
    mylikT1x <- t.y1[2] * log(t.pT1[1]) +  (t.n1[2] - t.y1[2]) * log(1 - t.pT1[1])
    mylikT1r <- t.y1[3] * log(t.pT1[2]) +  (t.n1[3] - t.y1[3]) * log(1 - t.pT1[2])
    
    mylik<-c(mylikL0, mylikL1, mylikT0, mylikT1x, mylikT1r)

    -sum(mylik)

  }


## fit model with pi in lik calc
fn.fitme4<-function(initpar, initdta)
  {

     p.pi0<-initdta[1]; t.n0<-initdta[2:3]; t.n1<-initdta[4:6]; t.y0<-initdta[7:8]; t.y1<-initdta[9:11]

     mydta<-c(NA, t.n0, t.n1, t.y0, t.y1); mypar<-initpar

     t.opt<-optim(par=mypar, fn=fn.lik4, mydta=mydta, method="BFGS")

##     print(t.opt)
     
     c(t.opt$value, exp(t.opt$par))
    
  }





###############


### simple approach to get first estimates
fn.initfit<-function(t.y0, t.y1, t.n0, t.n1)
  {
    p.theta1<-(t.y0[1] / t.y0[2]) * (t.n0[2] / t.n0[1])
    p.base2<-t.y1[3]/t.n1[3]
##    p.omega<-(t.y1[2]/t.n1[2]) / p.base2
    p.omega<-(t.y1[2]/t.n1[2]) / (p.base2*p.theta1)
    p.pi<-fn.solvepi(t.n0[2], t.y0[2], t.n0[2] - t.y0[2] - t.n1[2] -t.n1[3], p.omega, t.n1[2])
    p.pi1<-fn.fitu1(t.n0[1], t.y0[1], t.n0[1] - t.n1[1] - t.y0[1], p.omega, p.pi * t.n0[1]) / t.n1[1]
    p.cf<-p.pi1/ (t.n1[2] / sum(t.n1[2:3])) * (t.n1[1] / sum(t.n1[2:3]) ) ##correction factor
    p.theta2<-(t.y1[1] - p.cf * t.y1[2]) / (t.y1[3]) *  ( t.n1[3] / ( t.n1[1] * (1-p.pi1) ) )
  
    p.thetao<-(t.y0[1]+t.y1[1] - p.cf * t.y1[2]) / (t.y0[2]+ t.y1[3]) * (t.n0[2] / t.n0[1])
    p.itt<-(t.y0[1] + t.y1[1]) / (t.y0[2] + t.y1[2] + t.y1[3]) * (t.n0[2] / t.n0[1])

    myw<-c(t.n0[2], t.n1[3]) / (t.n0[2] + t.n1[3])
    myov<-c(p.theta1, p.theta2) %*% myw
    c(myov, p.theta1, p.theta2, p.omega, p.base2, p.pi, p.pi1, p.thetao, p.itt)
  }

## relative risk (ITT) profile likelihood cis fns
fn.lik.bin<-function(inpar, iny, inn)
  {
    log(inpar)*iny + log(1-inpar)*(inn-iny)
  }

fn.lik.bin2<-function(inpar1, intheta, iny, inn)
  {
     mylik<-fn.lik.bin(intheta*exp(inpar1), iny[1], inn[1]) +  fn.lik.bin(exp(inpar1), iny[2], inn[2])
     -mylik
  }

fn.fitrelrisk<-function(otheta, pstart, iny, inn)
  {
    t.op<-optim(pstart, fn.lik.bin2, intheta=otheta, iny=iny, inn=inn, method="BFGS")
    t.op$value
  }
fn.proflik<-function(otheta, pstart, iny, inn, ML, mychi)
  {
   thislik <- fn.fitrelrisk(otheta, pstart, iny, inn)
   thislik-ML -mychi
  }


#############################################################
## 2. Analysis
#############################################################

##DATA
##A. TRIAL
t.n0<-c(2463, 2459) ##n baseline L, T
t.n1<-c(2463-352-16-50,619, 2459-418-11-55-619) ##n period 1 -- updated using just number DFS + lost to fu

#B. OUTCOME
##B1. DFS
t.y0<-c(352,418) ## DFS L, T pre-switch
t.y1<-c(294, 58, 251) ## DFS post-switch L, T switch, T stay #12y


##FIT MODEL
p.initfit<-fn.initfit(t.y0, t.y1, t.n0, t.n1) 
p.theta<-p.initfit[1]
p.omega<-p.initfit[4]
p.base2<-p.initfit[5]
p.pi<-p.initfit[6]
p.base<-t.y0[2]/t.n0[2]

mypar<-c(p.base, p.base2, p.omega, p.theta)
myparmax<-1; myparscale1<- (log(mypar[1:2]) - log(myparmax))/(log(1-mypar[1:2]/myparmax))
myparmax<-1; myparscale2<- (log(mypar[3:4]) - log(myparmax))/(log(1-mypar[3:4]/myparmax))
myparinit<-c(myparscale1, myparscale2)

myparinit<-log(c(p.base, p.base2, p.omega, p.theta))
mydta<-c(p.pi, t.n0, t.n1, t.y0, t.y1)

myfit2<-fn.fitme2(myparinit, mydta)
myCI2<-fn.profCI2(myfit2, mydta, 0.05)
c(myfit2[5], myCI2)
p.initfit


##check constancy assumptions
myfit3<-fn.fitme3(c(myparinit, myparinit[4]), mydta)
-2*(myfit3[1]-myfit2[1])
postscript(file="fig_pl.eps", paper="special",width=6.0,height=6.0,horiz=F,pointsize=10)
par(mfrow=c(2,1))
myCI3a<-fn.profCI3(myfit3,  mydta,0.05,0, c(0.66,1.5),3.5)
title(main="(a)")
myCI3b<-fn.profCI3(myfit3,  mydta,0.05,1, c(0.66,1.5),3.5)
title(main="(b)")
dev.off()

fn.format<-function(ind, ndigit=2)
    {
                    format(round(ind,ndigit), nsmall=ndigit)
                }

fn.format(cbind(myfit3[5:6], rbind(myCI3a,myCI3b)))
-2*(myfit3[1]-myfit2[1])
1-pchisq(-2*(myfit3[1]-myfit2[1]) ,1)


myfit4<-fn.fitme4(c(myparinit, myparinit[1]), mydta)
-2*(myfit4[1]-myfit2[1])

##myfit5<-fn.fitme5(c(myparinit, myparinit[1], myparinit[4]), mydta)
##-2*(myfit5[1]-myfit3[1])

p.pi<-fn.solvepi(t.n0[2], t.y0[2], t.n0[2] - t.y0[2] - t.n1[2] -t.n1[3], myfit2[4], t.n1[2])

p.pi1<-fn.fitu1(t.n0[1], t.y0[1], t.n0[1] - t.n1[1] - t.y0[1], myfit2[4], p.pi * t.n0[1]) / t.n1[1]

##itt ci.
myarm<-c(rep(1,t.n0[1]), rep(0, t.n0[2]))

myout1<-c(rep(1, t.y0[1]+t.y1[1]), rep(0, t.n0[1]-t.y0[1]-t.y1[1]))

myout0<-c(rep(1, t.y0[2]+t.y1[2]+t.y1[3]), rep(0, t.n0[2]-t.y0[2]-t.y1[2]-t.y1[3]))

myy<-c(myout1, myout0)

##odds ratio
t.glm<-glm(myy~myarm, family=binomial(link="logit"))
exp(confint(t.glm))
exp(coef(t.glm))

##rel risk
iny<-c(t.y0[1]+t.y1[1], t.y0[2]+t.y1[2]+t.y1[3]);
inn<-t.n0
pstart<-log(iny[2]/inn[2])
otheta<- (iny[1]/inn[1]) / exp(pstart)

myML<-fn.fitrelrisk(otheta, pstart, iny, inn)
myalpha<-0.05;mychi<-qchisq(1-myalpha, 1) / 2
fn.proflik(otheta+0.1, pstart, iny, inn, myML, mychi)

myCIU<-uniroot(fn.proflik, interval=c(0,0.5)+otheta, pstart=pstart, mychi=mychi, ML=myML, inn=inn, iny=iny)
myCIL<-uniroot(fn.proflik, interval=c(-0.5,0)+otheta, pstart=pstart, mychi=mychi, ML=myML, inn=inn, iny=iny)
c(otheta, myCIL$root, myCIU$root)



