################################################################################################w
## R functions for selective xover estimation
## Adam Brentnall
## This has been worked on since the code used in the paper. 
## 12th April 2018
################################################################################################w

##MAIN ROUTINES
## Fit and estimate CIs, stratified model
fn.main.stratified<-function(infile, doci=TRUE)
  {

##    infile<-"simin/debug.csv" ##debug
    ##source("fitfns.r");    dyn.load("Rmain.so"); library("boot"); library("survival")
    
      ##1. LOAD DATA
    mydta<-  fn.load.data(infile)
  ##  mydta2<-  fn.load.data.breakties(infile)   

    ##2. FIT
##    myfit2<-fn.main.stratified.fit(mydta2)
    myfit<-fn.main.stratified.fit(mydta)


    ##3. CIs
    if(doci==TRUE)
      {
        mycis<-fn.main.stratified.ci(mydta, myfit[[2]], myfit[[3]], inR=99)
        myout<-    list(itt=myfit[[1]], plest=myfit[[2]][1], mycis)
      }
    else
      {
        myout<-list(myfit[[1]], myfit[[2]][1], myfit[[3]][1])
##        myout<-list(myfit[[1]], myfit[[2]][1], myfit[[3]][1], myfit2[[1]], myfit2[[2]][1], myfit2[[3]][1])
      }
    myout
  }


##Fit only
fn.main.stratified.fit<-function(mydta)
  {
    ##Input: data 
    
    ##1. INITIALISE starting parameter fit
    mypstart<-fn.initialise(mydta)

    #################################################################
    ##2. FIT STRATIFIED
    myopt.st.pl<- fn.PL.fit.dbg(mypstart, mydta)
    myopt.st<-optim(myopt.st.pl[[1]], fn=fn.FULL.dbg,  method="BFGS", ind=mydta) ##use PL as start
    myopt.st.obj<- fn.FULL.all.dbg(myopt.st$par, mydta)    

    list(itt=mypstart[2], myopt.st.pl, myopt.st.obj)
  }

#CI's.
fn.main.stratified.ci<-function(mydta, myopt.st.pl, myopt.st.obj, inR=49)
  {
##INPUT: data, partial likelihood fitted object, full likelihood fitted object, number of bootstraps for PL

    ##3. CIs
    ##3A. PL bootstrap CI
    mybsci<- fn.plboot(mydta, myopt.st.pl[[1]], inR)
    outpc<-sapply(1:2, function(idx) boot.ci(mybsci, index=idx, type="perc")$percent[4:5])

    ##3B. Profile likelihood CI
    mylow<-myopt.st.obj[[1]][2]- 3*(myopt.st.pl[[1]][2] - outpc[1,2]) ##look three times either side of pl interval
    myhi<-myopt.st.obj[[1]][2]+ 3*(-myopt.st.pl[[1]][2] + outpc[2,2])
    outproflikci<-fn.FULL.proflik(myopt.st.obj[[1]], myopt.st.obj[[2]], mydta, inalp=0.05, inlow=mylow, inup=myhi)

    list(ph=outpc, fl=outproflikci)
  }



##Use this if know that no strata
fn.main.xstratified<-function(infile)
  {

    ##source("fitfns.r");    dyn.load("Rmain.so"); library("boot"); library("survival")

    ##1. LOAD DATA
    ##dedug -- gendata
    infile<-"simin/sim1.csv";     insp<-fn.simscen(2000);     fn.gensamp(1,insp[[2]],infile)
    mydta<-  fn.load.data(infile)
##    mydta<-  fn.load.data("demoinxs.csv")
    
  
    
    ##2. INITIALISE starting parameter fit
    mypstart<-fn.initialise(mydta)

    #################################################################
    ##3. FIT UN STRATIFIED
    myopt.pl<-fn.PL.fit(mypstart, mydta )
    myopt<-optim(myopt.pl[[1]], fn=fn.FULL,  method="BFGS", ind=mydta)
    myopt.obj<- fn.FULL.all(myopt$par, mydta)    

    ##3. CIs
    ##3A. PL bootstrap CI
    inR<-9
    mybsci<- fn.plboot(mydta, myopt.st.pl[[1]], inR, coxstrata=FALSE)
    outpc<-sapply(1:2, function(idx) boot.ci(mybsci, index=idx, type="perc")$percent[4:5])

    ##3B. Profile likelihood CI
    mylow<-myopt.obj[[1]][2]- 3*(myopt.pl[[1]][2] - outpc[1,2]) ##look three times either side of pl interval
    myhi<-myopt.obj[[1]][2]+ 3*(-myopt.pl[[1]][2] + outpc[2,2])
    outproflikci<-fn.FULL.proflik(myopt.obj[[1]], myopt.obj[[2]], mydta, inalp=0.05, inlow=mylow, inup=myhi, coxstrata=FALSE)

    list(itt=mypstart[2], ph=outpc, fl=outproflikci)
    
  }



##########################
##Load data
##Loads csv file.
##FORMAT:
## ROWS - number of participants plus one row (header)
## COLUMNS - 7 columns:
##   1 id;
##   2 time to event;
##   3 event indicator (1=event, 0 = censored);
##   4 xover group (0 - unknown, 1 = don't xover, 2 = xover);
##   5 arm (1 or 2);
##   6 potential xover time (calculatd for everyone = calender time at unblinding - entry);
##   7 strata (code as 1,2,3,...,number of baseline hazard strata in Cox model).
###########################
fn.load.data<-function(inpath)
  {
    ##INPUT
    ##inpath: full path to csv file

    ##debug
##    inpath<-"demoin.csv"
    
    ##load data
    myin<-read.csv(inpath, header=TRUE)

    ##order
    myorder<-order(myin[,7], myin[,2])
    myin<-myin[myorder,]                 
    
    ##derive other fields needed for program
    myout<-list()

    ##number of participants
    mynT<-nrow(myin)
    myn1<-sum(myin[,5]==1)
    myn2<-sum(myin[,5]==2)

    ##covariates
    mym<-1 ##number of x's (dummy)
    mydatax<-rep(0,mynT)
    
    ##strata
    nstrata<-max(myin[,7])
    mystrata<-cumsum(c(0,hist(myin[,7],plot=F, breaks=seq(0.5, nstrata+0.5, by=1))$count)) #positions of strata

    ## event time point before than individual switched, by strata
    datasidx<-unlist(lapply(1:nstrata, function(idx) sapply(seq(mystrata[idx]+1, mystrata[idx+1]) , function(idy) pmin(diff(mystrata[idx:(idx+1)])-1,sum(myin[seq(mystrata[idx]+1, mystrata[idx+1]),2]< myin[idy,6])))))
    
    ##put together
    myall<-list(mynT=mynT,
                mydatat=myin[,2],
                mydatad=myin[,3],
                mydatag=myin[,4],
                mydatar=myin[,5],
                mydatas=myin[,6],
                mydatasidx=datasidx,
                mystrata=mystrata,
                mynstrata=nstrata,
                
                mym=1, ##number covariates
                mydatax=rep(0,mynT), #covariates (none used here)
                mypl=0, #partial likelihood,
                mylik=0, #likliehood                
                myp1=rep(0,myn1), #other terms returned from fitting - proportion of switcher through time estimates
                myp3=rep(0,myn2),
                mypitime=rep(0,mynT))

    ##return
    myall
               
  }

##break ties by adding a small number. Needed for baselinehazard full likelihood
fn.load.data.breakties<-function(inpath)
  {
    ##INPUT
    ##inpath: full path to csv file

    ##debug
##    inpath<-"demoin.csv"
    
    ##load data
    myin<-read.csv(inpath, header=TRUE)

    ##order
    counter<-0; mySMALLNUM<-min(myin[,2]/100000, 0.00000000001)
    mytiet<-myin[,2]
    for(idx in 2:nrow(myin))
      {

         if(myin[idx,2] == myin[idx-1,2])
           {
             counter<-counter+1
             mytiet[idx]<-myin[idx-1,2]+mySMALLNUM*counter
           }
         else
           {            
             counter<-0
           }
       }

    myin[,2]<-mytiet
    
    myorder<-order(myin[,7], myin[,2])
    myin<-myin[myorder,]                 
    
    ##derive other fields needed for program
    myout<-list()

    ##number of participants
    mynT<-nrow(myin)
    myn1<-sum(myin[,5]==1)
    myn2<-sum(myin[,5]==2)

    ##covariates
    mym<-1 ##number of x's (dummy)
    mydatax<-rep(0,mynT)
    
    ##strata
    nstrata<-max(myin[,7])
    mystrata<-cumsum(c(0,hist(myin[,7],plot=F, breaks=seq(0.5, nstrata+0.5, by=1))$count)) #positions of strata

    ## event time point before than individual switched, by strata
    datasidx<-unlist(lapply(1:nstrata, function(idx) sapply(seq(mystrata[idx]+1, mystrata[idx+1]) , function(idy) pmin(diff(mystrata[idx:(idx+1)])-1,sum(myin[seq(mystrata[idx]+1, mystrata[idx+1]),2]< myin[idy,6])))))

    
    ##put together
    myall<-list(mynT=mynT,
                mydatat=myin[,2],
                mydatad=myin[,3],
                mydatag=myin[,4],
                mydatar=myin[,5],
                mydatas=myin[,6],
                mydatasidx=datasidx,
                mystrata=mystrata,
                mynstrata=nstrata,
                
                mym=1, ##number covariates
                mydatax=rep(0,mynT), #covariates (none used here)
                mypl=0, #partial likelihood,
                mylik=0, #likliehood                
                myp1=rep(0,myn1), #other terms returned from fitting - proportion of switcher through time estimates
                myp3=rep(0,myn2),
                mypitime=rep(0,mynT))

    ##return
    myall
               
  }


##Sensible starting estimates. Cox ITT for trt effect.
fn.initialise<-function(sim.dta)
  {
        mystrata<-rep(1:(sim.dta$mynstrata), diff(sim.dta$mystrata))
        
        t.fit<-coxph(Surv(sim.dta$mydatat, sim.dta$mydatad) ~ sim.dta$mydatar+strata(mystrata))
        
        t.inc<-sim.dta$mydatar==1
        
        t.fit2<-coxph(Surv(sim.dta$mydatat[t.inc], sim.dta$mydatad[t.inc]) ~ (sim.dta$mydatag[t.inc]==2) + strata(mystrata[t.inc]))
        
        myout<-c(coef(t.fit2)/coef(t.fit), coef(t.fit))
        
        myout
  }

## Bootstrap confidence intervals from partial likelihood
fn.plboot<-function(sim.dta, myfit, inR=99, coxstrata=TRUE)
  {

    ##PL fits from stratified model
    fn.plboot.stat<-function(myrids, inidx, sim.dta, myfit)
      {
        tidx<-myrids[inidx]       
       
        ## create bootstrap data
        mydta<-sim.dta

        #time
        mydta$mydatat<-sim.dta$mydatat[tidx];

        ##order
        myord.l<-sapply(1:(sim.dta$mynstrata), function(idx) sim.dta$mystrata[idx]+ order(mydta$mydatat[seq(sim.dta$mystrata[idx]+1, sim.dta$mystrata[idx+1])]))       
        myord<-unlist(myord.l)
              
        #order everything appropriately
        mydta$mydatat<-mydta$mydatat[myord]
        mydta$mydatad<-sim.dta$mydatad[tidx[myord]]
        mydta$mydatax<-sim.dta$mydatax[tidx[myord]]; mydta$mydatar<-sim.dta$mydatar[tidx[myord]]
        mydta$mydatas<-sim.dta$mydatas[tidx[myord]]; mydta$mydatag<-sim.dta$mydatag[tidx[myord]]

        if(coxstrata)
          {
            temp.pl.boot<- fn.PL.fit.dbg(myfit, mydta)
          }
        else
          {
             temp.pl.boot<- fn.PL.fit(myfit, mydta)
          }

       temp.pl.boot[[1]]

      }

##debug
##    sim.dta<-mydta; myfit<-myopt.st.pl[[1]]; inR<-2
    
    myrids<-c(seq(1,sim.dta$mynT))
    mystrata<-rep(1:(sim.dta$mynstrata), diff(sim.dta$mystrata))
    mybstrata<-  sim.dta$mynstrata*(sim.dta$mydatar-1) + mystrata

    t.boot<-boot(data=myrids, statistic=fn.plboot.stat, strata=mybstrata, sim.dta=sim.dta, myfit=myfit, R=inR)

    t.boot

  }




## Used by profile likelihood. Treatment effect fixed, pass in ind as $trteff
fn.FULL.trtfixed<-function(inpar, ind)
      {

        fitpar<-c(inpar[1], ind$trteff) ##no covariates yet
            
        temp<-.C("testerfull_strata", as.double(fitpar), as.double(ind$mylik), as.double(ind$mydatat), as.integer(ind$mydatad), as.double(ind$mydatax), as.integer(ind$mydatag), as.integer(ind$mydatar), as.integer(ind$mynT), as.integer(ind$mym), as.double(ind$mydatas), as.double(ind$myp1), as.double(ind$mypitime), as.double(ind$mypl), as.integer(ind$mydatasidx), as.double(ind$myp3), as.integer(ind$mystrata), as.integer(ind$mynstrata))
        print(paste("FULL LIK",round(temp[[2]],2), "Prof par:",round(exp(fitpar[1]),4) ))
       temp[[2]]
      }

## Used by profile likelihood. Treatment effect fixed, pass in ind as $trteff
## Not stratified baseline
fn.FULLxs.trtfixed<-function(inpar, ind)
      {

        fitpar<-c(inpar[1], ind$trteff) ##no covariates yet
            
        temp<-.C("testerfull", as.double(fitpar), as.double(ind$mylik), as.double(ind$mydatat), as.integer(ind$mydatad), as.double(ind$mydatax), as.integer(ind$mydatag), as.integer(ind$mydatar), as.integer(ind$mynT), as.integer(ind$mym), as.double(ind$mydatas), as.double(ind$myp1), as.double(ind$mypitime), as.double(ind$mypl), as.integer(ind$mydatasidx), as.double(ind$myp3), as.integer(ind$mystrata), as.integer(ind$mynstrata))
        print(paste("FULL LIK",round(temp[[2]],2), "profile par:",exp(fitpar[1])))
       temp[[2]]
      }


## profile likelihood confidence interval
## mlpar is MLE all pars, mllik is corresponding log likelihood, inalp is for confint
## likelihood fit, IN: parameters, data; OUT: fit 
fn.FULL.proflik<-function(mlpar, mllik, ind, inalp=0.05, inlow=log(0.1), inup=log(10), coxstrata=TRUE)
      {
        ind$trteff<-mlpar[2]
        mychi<-qchisq(1-inalp, 1) 

        fn.rootfind.1<-function(intrt, inother, mllik, mychi, ind, coxstrata)
          {
            inpar<-c(inother, intrt)
            ind$trteff<-intrt
            if(coxstrata==TRUE)
              {
                thisfit<-optim(inpar, fn=fn.FULL.trtfixed,  method="BFGS", ind=ind)
              }
            else
              {
                thisfit<-optim(inpar, fn=fn.FULLxs.trtfixed,  method="BFGS", ind=ind)
              }
            thislik<-thisfit[[2]]
            myout<-(2*(mllik - thislik) + mychi)
            myout
          }

##       mylower<-log(0.1) #lowest possible -- to investigate
##        myupper<-log(10) #lowest possible -- to investigate changing 
      
       ci.lower<-uniroot(fn.rootfind.1, interval=c(inlow, mlpar[2]) , inother=mlpar[1], mllik=mllik, mychi=mychi, ind=ind, coxstrata=coxstrata)
       ci.upper<-uniroot(fn.rootfind.1, interval=c(mlpar[2], inup) , inother=mlpar[1], mllik=mllik, mychi=mychi, ind=ind, coxstrata=coxstrata)

       myout<-(c(ci.lower[[1]], mlpar[2],ci.upper[[1]]))
       myout
      }


####################################################
## LIKELIHOOD FUNCTIONS AND FITTING VIA .C CALLS ###
####################################################

##1. STRATIFIED LIKELIHOODS

##partial likelihood, all the object
fn.PL.ALL.dbg<-function(inpar, ind)
      {

        temp<-.C("rint_pl_strata", as.double(inpar), as.double(ind$mylik), as.double(ind$mydatat), as.integer(ind$mydatad), as.double(ind$mydatax), as.integer(ind$mydatag), as.integer(ind$mydatar), as.integer(ind$mynT), as.integer(ind$mym), as.double(ind$mydatas), as.double(ind$myp1), as.double(ind$mypitime), as.double(ind$mypl), as.integer(ind$mystrata), as.integer(ind$mynstrata))
        temp
      }

fn.PL.fit.dbg<-function(inpar, ind)
      {

        temp<-.C("testerpl_strata", as.double(inpar), as.double(ind$mylik), as.double(ind$mydatat), as.integer(ind$mydatad), as.double(ind$mydatax), as.integer(ind$mydatag), as.integer(ind$mydatar), as.integer(ind$mynT), as.integer(ind$mym), as.double(ind$mydatas), as.double(ind$myp1), as.double(ind$mypitime), as.double(ind$mypl), as.integer(ind$mystrata), as.integer(ind$mynstrata))
        temp
      }

## Full likelihood fit, IN: parameters, data; OUT: fit 
fn.FULL.dbg<-function(inpar, ind)
      {
                   
            temp<-.C("testerfull_strata", as.double(inpar), as.double(ind$mylik), as.double(ind$mydatat), as.integer(ind$mydatad), as.double(ind$mydatax), as.integer(ind$mydatag), as.integer(ind$mydatar), as.integer(ind$mynT), as.integer(ind$mym), as.double(ind$mydatas), as.double(ind$myp1), as.double(ind$mypitime), as.double(ind$mypl), as.integer(ind$mydatasidx), as.double(ind$myp3), as.integer(ind$mystrata), as.integer(ind$mynstrata))
            print(paste("FULL LIK",round(temp[[2]],2), "pars:",exp(inpar[1]), exp(inpar[2])))
            temp[[2]]
      }

fn.FULL.all.dbg<-function(inpar, ind)
      {

#        thispar<-log(exp(inpar)/(1+exp(inpar)) * 100)
        temp<-.C("testerfull_strata", as.double(inpar), as.double(ind$mylik), as.double(ind$mydatat), as.integer(ind$mydatad), as.double(ind$mydatax), as.integer(ind$mydatag), as.integer(ind$mydatar), as.integer(ind$mynT), as.integer(ind$mym), as.double(ind$mydatas), as.double(ind$myp1), as.double(ind$mypitime), as.double(ind$mypl), as.integer(ind$mydatasidx), as.double(ind$myp3), as.integer(ind$mystrata), as.integer(ind$mynstrata))
        temp
      }


#########################
### 2. Not stratified ###
#########################

##partial likelihood fit, IN: parameters, data; OUT: fit 
fn.PL<-function(inpar, ind)
      {

        temp<-.C("rint_pl", as.double(inpar), as.double(ind$mylik), as.double(ind$mydatat), as.integer(ind$mydatad), as.double(ind$mydatax), as.integer(ind$mydatag), as.integer(ind$mydatar), as.integer(ind$mynT), as.integer(ind$mym), as.double(ind$mydatas), as.double(ind$myp1), as.double(ind$mypitime), as.double(ind$mypl))
        print(paste("PL",round(temp[[13]],2), "pars:", exp(inpar[1]), exp(inpar[2])))
        temp[[13]]
      }
#partial likelihood, all the object
fn.PL.ALL<-function(inpar, ind)
      {

        temp<-.C("rint_pl", as.double(inpar), as.double(ind$mylik), as.double(ind$mydatat), as.integer(ind$mydatad), as.double(ind$mydatax), as.integer(ind$mydatag), as.integer(ind$mydatar), as.integer(ind$mynT), as.integer(ind$mym), as.double(ind$mydatas), as.double(ind$myp1), as.double(ind$mypitime), as.double(ind$mypl))
        temp
      }


##partial likelihood fit, IN: parameters, data; OUT: fit 
fn.PL.fit<-function(inpar, ind)
      {

        temp<-.C("testerpl", as.double(inpar), as.double(ind$mylik), as.double(ind$mydatat), as.integer(ind$mydatad), as.double(ind$mydatax), as.integer(ind$mydatag), as.integer(ind$mydatar), as.integer(ind$mynT), as.integer(ind$mym), as.double(ind$mydatas), as.double(ind$myp1), as.double(ind$mypitime), as.double(ind$mypl))
        temp
      }


## FULL LIKELIHOOD not stratified
fn.FULL<-function(inpar, ind)
      {
                   
            temp<-.C("testerfull", as.double(inpar), as.double(ind$mylik), as.double(ind$mydatat), as.integer(ind$mydatad), as.double(ind$mydatax), as.integer(ind$mydatag), as.integer(ind$mydatar), as.integer(ind$mynT), as.integer(ind$mym), as.double(ind$mydatas), as.double(ind$myp1), as.double(ind$mypitime), as.double(ind$mypl), as.integer(ind$mydatasidx), as.double(ind$myp3))
            print(paste("FULL LIK",round(temp[[2]],2), "pars:",exp(inpar[1]), exp(inpar[2])))
            temp[[2]]
      }

fn.FULL.all<-function(inpar, ind)
      {

        temp<-.C("testerfull", as.double(inpar), as.double(ind$mylik), as.double(ind$mydatat), as.integer(ind$mydatad), as.double(ind$mydatax), as.integer(ind$mydatag), as.integer(ind$mydatar), as.integer(ind$mynT), as.integer(ind$mym), as.double(ind$mydatas), as.double(ind$myp1), as.double(ind$mypitime), as.double(ind$mypl), as.integer(ind$mydatasidx), as.double(ind$myp3))
        temp
      }


####################################
####################################
##  BELOW: simulation functions  ###
####################################
####################################
fn.genbigsim<-function(inseed = 101, infile="simin/debug.txt")
  {

    ##debug
   ##inseed<-1032; infile="simin/debug.txt"
    
    
    my2par<-c(1.5, 0.7, 0.05, 0.1) #tau, theta, cen before, cen after

    ##strata 1: 4 arm, no chemo
    myn<-1160*2
    myst1<-list(myn=myn, mypi=0.5, mypar=list(tau=my2par[1], theta=my2par[2], base=fn.base(c(1,my2par[2:4]))), mycen=fn.cen(c(1,my2par[2:4])), lsprob=0.5, lstime=9999999999)

    ##strata 2: 4 arm, chemo
    myst2<-myst1
    myst2$base=myst1$base*1.1
    myst2$myn<-390*2

    ##strata 3: 2 arm, no chem0
    myst3<-myst1
    myst3$base=myst1$base*0.8
    myst3$myn<-700*2

    ##strata 4: 2 arm, chemo
    myst4<-myst1
    myst4$base=myst1$base*0.9
    myst4$myn<-200*2

    insp.l<-list(myst1, myst2, myst3, myst4)

    myout<-fn.strat.gensamp2(inseed, insp.l, infile)  
    
  }

### simulation scenarios
fn.simscen<-function(myn)
  {
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

    insp
}


##used to obtain start end points of sim for censoring
fn.base<-function(ind)
  {
    intau<-ind[1]; intheta<-ind[2]; incen<-ind[3]
    tout<--log(1-incen)/(intau*intheta)
  }
fn.cen<-function(ind)
  {
    intau<-ind[1]; intheta<-ind[2]; incen1<-ind[3]; incen2<-ind[4]
    myrate<-intau * intheta * fn.base(ind)
    tout<-qexp( min(0.9999999,incen1*2+incen2), myrate)
  }



##write simulated data to output for input to main prog
fn.writecsvinout<-function(ind, infile)
  {
##    ind<-sim.dta ##debug
    if(is.null(ind$mynstrata)==FALSE      ){
        myst<-rep(seq(1,ind$mynstrata), diff(ind$mystrata))
      }    else      {
        myst<-rep(1,ind$mynT)
      }

     myout<-data.frame(id=seq(1,ind$mynT), 
     t=ind$mydatat,
     d=ind$mydatad,
     g=ind$mydatag,
     r=ind$mydatar,
     sx=ind$mydatas,
     st=myst)

   write.csv(myout, file=infile, row.names=FALSE)

  }

## generate sample not stratified
fn.gensamp<-function(inseed, insp, infile)
  {
    print(paste("RSEED",inseed))
    set.seed(inseed)   
    sim.dta.all<-fn.simdata(insp)
    sim.dta<-sim.dta.all[[1]]
    fn.writecsvinout(sim.dta, infile=infile)
    
  }


## generate stratified sample with two strata, baseline haz set in function
fn.strat.gensamp<-function(inseed, insp, infile)
  {
   print(paste("RSEED",inseed))
    set.seed(inseed)
## insp<-insp1[[8]] ##debug
   
   sim.dta.all1<-fn.simdata.dbg(insp)
   ##add strat1a
   insp$mypar$base<-insp$mypar$base*1.2
   sim.dta.all2<-fn.simdata.dbg(insp)

   sim.dta.all1$obs$mydatasidx<-sapply(1:sim.dta.all1$obs$mynT, function(idx) pmin(sim.dta.all1$obs$mynT-1,sum(sim.dta.all1$obs$mydatat<sim.dta.all1$obs$mydatas[idx])))
   sim.dta.all2$obs$mydatasidx<-sapply(1:sim.dta.all2$obs$mynT, function(idx) pmin(sim.dta.all2$obs$mynT-1,sum(sim.dta.all2$obs$mydatat<sim.dta.all2$obs$mydatas[idx])))

   sim.dta.all<-list(obs = list(mylik=as.double(-99), mydatat=c(sim.dta.all1$obs$mydatat, sim.dta.all2$obs$mydatat),
                       mydatad=c(sim.dta.all1$obs$mydatad, sim.dta.all2$obs$mydatad),
                       mydatax=c(sim.dta.all1$obs$mydatax, sim.dta.all2$obs$mydatax),
                       mydatag=c(sim.dta.all1$obs$mydatag, sim.dta.all2$obs$mydatag),
                       mydatar=c(sim.dta.all1$obs$mydatar, sim.dta.all2$obs$mydatar),
                       mynT = sim.dta.all1$obs$mynT + sim.dta.all2$obs$mynT,
                       mym=sim.dta.all1$obs$mym,
                       myp1=c(sim.dta.all1$obs$myp1, sim.dta.all2$obs$myp1),
                       mypitime=c(sim.dta.all1$obs$mypitime, sim.dta.all2$obs$mypitime),
                       myp3=c(sim.dta.all1$obs$myp3, sim.dta.all2$obs$myp3),
                       mypl=0,
                       mydatas=c(sim.dta.all1$obs$mydatas, sim.dta.all2$obs$mydatas),
                       mydatasidx= c(sim.dta.all1$obs$mydatasidx, sim.dta.all2$obs$mydatasidx)),
                     unobs=list(c(unlist(sim.dta.all1$unobs), unlist(sim.dta.all2$unobs))))

   myord<-order(sim.dta.all$obs$mydatat)
   sim.dta.all.sort<-list(obs = list(mylik=as.double(-99), mydatat=c(sim.dta.all1$obs$mydatat, sim.dta.all2$obs$mydatat)[myord],
                       mydatad=c(sim.dta.all1$obs$mydatad, sim.dta.all2$obs$mydatad)[myord],
                       mydatax=c(sim.dta.all1$obs$mydatax, sim.dta.all2$obs$mydatax)[myord],
                       mydatag=c(sim.dta.all1$obs$mydatag, sim.dta.all2$obs$mydatag)[myord],
                       mydatar=c(sim.dta.all1$obs$mydatar, sim.dta.all2$obs$mydatar)[myord],
                       mynT = sim.dta.all1$obs$mynT + sim.dta.all2$obs$mynT,
                       mym=sim.dta.all1$obs$mym,
                       myp1=c(sim.dta.all1$obs$myp1, sim.dta.all2$obs$myp1),
                       mypitime=c(sim.dta.all1$obs$mypitime, sim.dta.all2$obs$mypitime),
                       myp3=c(sim.dta.all1$obs$myp3, sim.dta.all2$obs$myp3),
                       mypl=0,
                       mydatas=c(sim.dta.all1$obs$mydatas, sim.dta.all2$obs$mydatas)[myord],
                       mydatasidx= c(sim.dta.all1$obs$mydatasidx, sim.dta.all2$obs$mydatasidx)[myord]),                            
                     unobs=list(c(unlist(sim.dta.all1$unobs), unlist(sim.dta.all2$unobs))[myord] ))

   
   sim.dta<-sim.dta.all.sort[[1]]
   sim.z<-sim.dta.all.sort[[2]] 

   sim.dta$mystrata<-c(0,sim.dta$myn/2, sim.dta$myn)
   sim.dta$mynstrata<-2
   
    fn.writecsvinout(sim.dta, infile)    
  }



## generate and save stratified sample for length(inspl) strata
fn.strat.gensamp2<-function(inseed, insp.l, infile)
  {
    ##debug
#    my2par<-c(1.5, 0.7, 0.05, 0.1) #tau, theta, cen before, cen after
#    myn<-1000
#    myst1<-list(myn=myn, mypi=0.5, mypar=list(tau=my2par[1], theta=my2par[2], base=fn.base(c(1,my2par[2:4]))), mycen=fn.cen(c(1,my2par[2:4])), lsprob=0.5, lstime=9999999999)
#    myst2<-myst1
#    myst2$base=myst1$base*1.1
#    myst2$myn<-500
 #   insp.l<-list(myst1, myst2)
    ## end debug
    
   print(paste("RSEED",inseed))
    set.seed(inseed)
   
##    sim.dta.all.l<-lapply(insp.l, fn.simdata.bigstrata) ##no ties
    sim.dta.all.l<-lapply(insp.l, function(ind) fn.simdata.bigstrata.round(ind, 2)) ##rounding error

    myinc<-c(2,3,5,6,13)
    mystat<-    sapply(myinc, function(idx) unlist(lapply(sim.dta.all.l, function(ind) ind[[1]][[idx]])))
    mystrata<-rep(seq(1, length(sim.dta.all.l)), unlist(lapply(sim.dta.all.l, function(ind) ind[[1]]$mynT)) )
    myorder<-order(mystat[,1])

    mydta<-cbind(mystat, mystrata)[myorder,]

    myout<-data.frame(id=seq(1,nrow(mydta)),
     t=pmax(0.000001,mydta[,1]), 
     d=mydta[,2],
     g=mydta[,3],
     r=mydta[,4],
     sx=pmax(0.000001,mydta[,5]),
     st=mydta[,6])

   write.csv(myout, file=infile, row.names=FALSE)

    myout

   
  }

##used to generate stratified samples, where some have large switch time
fn.simdata.bigstrata<-function(insp)
  {
## debug
    ##    myn<-1000;  mypi<-0.25; mypar<-list(tau=0.5, theta=0.7, base=1)    
##    insp<-list(myn=1000, mypi=0.25, mypar=list(tau=0.5, theta=0.7, base=1), mycen=10)

    ## possible to change/add  later on
    tbin<-1; tbin2<-2;   mym<-1 #tbin2 = tbin then all switch at same tbin time

    
    ##censoring
#    mymaxt1<-insp$mycen*runif(insp$myn);  mymaxt2<-insp$mycen*runif(insp$myn)
    mymaxt1<-insp$mycen; mymaxt2<-insp$mycen

    ##switcher?
    myz1<-runif(insp$myn)<insp$mypi
    mys1<-tbin + runif(insp$myn)*(tbin2-tbin) #switch in second period ##random switching

    ##some with large switch time
    mylarges<-runif(insp$myn)>insp$lsprob
    mys1[mylarges==TRUE]<-insp$lstime #large switch time
    
    ##arm 1
    myt1as<-rexp(insp$myn, rate=( ( (1-myz1) + myz1*(insp$mypar$theta*insp$mypar$tau)) * insp$mypar$base))
    myt1ps<-rexp(insp$myn, rate=( ( (1-myz1) + myz1*insp$mypar$tau) * insp$mypar$base));  myt1ps<-pmin(myt1ps, mys1)
    myt1<-myt1ps;     myt1[myt1ps==mys1]<-(myt1ps+myt1as)[myt1ps==mys1]
    myd1<-myt1<mymaxt1 ;  myt1x<-pmin(myt1,mymaxt1)    ##censoring

    ##arm2 
##    myz2<-myz1;  mys2<-mys1     #switchers identical number to arm 1
    myz2<-runif(insp$myn)<insp$mypi
    mys2<-tbin + runif(insp$myn)*(tbin2-tbin) #switch in second period ##random switching

    mylarges2<-runif(insp$myn)>insp$lsprob
    mys2[mylarges2==TRUE]<-insp$lstime #large switch time
    
    myt2<-rexp(insp$myn, rate=((insp$mypar$theta *  (1-myz2+myz2*insp$mypar$tau) * insp$mypar$base))) #failure time, arm 2 
    myd2<-myt2<mymaxt2;    myt2x<-pmin(myt2,mymaxt2) ##censoring

    mydatad<-as.integer(c(myd1,myd2))
    mydatat<-c(myt1x, myt2x); mydatas<-c(mys1, mys2);  mydatar<-rep(1:2,each=insp$myn); mydatax<-rep(1,insp$myn*2)
    mydatag<-c((1+myz1)*(myt1x>=mys1), rep(0,insp$myn)) #0 = unknown, 1=don't switch, 2=switch
   
#reorder everything by time
    myord<-order(mydatat)
    mydatat<-mydatat[myord];  mydatas<-mydatas[myord]
    mydatax<-mydatax[myord];  mydatad<-mydatad[myord]
    mydatag<-mydatag[myord];  mydatar<-mydatar[myord]
    mynT<-2*insp$myn

    mydataz<-c(myz1, myz2)[myord]   
    
    myfitobj<-list(obs = list(mylik=as.double(-99), mydatat=mydatat, mydatad=mydatad, mydatax=mydatax, mydatag=mydatag, mydatar=mydatar, mynT=mynT, mym=mym, myp1=rep(0,insp$myn), mypitime=rep(0,insp$myn), myp3=rep(0,insp$myn), mypl=0, mydatas=mydatas), unobs=list(mydataz))

    myfitobj


  }


##used to generate stratified samples, where some have large switch time
## rounded to x dp
fn.simdata.bigstrata.round<-function(insp, mypre=2)
  {
## debug
    ##    myn<-1000;  mypi<-0.25; mypar<-list(tau=0.5, theta=0.7, base=1)    
##    insp<-list(myn=1000, mypi=0.25, mypar=list(tau=0.5, theta=0.7, base=1), mycen=10)

    ## possible to change/add  later on
    tbin<-1; tbin2<-2;   mym<-1 #tbin2 = tbin then all switch at same tbin time
    
    ##censoring
#    mymaxt1<-insp$mycen*runif(insp$myn);  mymaxt2<-insp$mycen*runif(insp$myn)
    mymaxt1<-round(insp$mycen, mypre); mymaxt2<-round(insp$mycen, mypre)

    ##switcher?
    myz1<-runif(insp$myn)<insp$mypi
    mys1<-tbin + runif(insp$myn)*(tbin2-tbin) #switch in second period ##random switching
    mys1<-round(mys1, mypre) ##rounding
    
    ##some with large switch time
    mylarges<-runif(insp$myn)>insp$lsprob
    mys1[mylarges==TRUE]<-insp$lstime #large switch time
    
    ##arm 1
    myt1as<-rexp(insp$myn, rate=( ( (1-myz1) + myz1*(insp$mypar$theta*insp$mypar$tau)) * insp$mypar$base))
    myt1ps<-rexp(insp$myn, rate=( ( (1-myz1) + myz1*insp$mypar$tau) * insp$mypar$base));  myt1ps<-pmin(myt1ps, mys1)
    myt1<-myt1ps;     myt1[myt1ps==mys1]<-(myt1ps+myt1as)[myt1ps==mys1]
    myt1<-round(myt1,mypre) ##rounding
    myd1<-myt1<mymaxt1 ;  myt1x<-pmin(myt1,mymaxt1)    ##censoring

    ##arm2 
##    myz2<-myz1;  mys2<-mys1     #switchers identical number to arm 1
    myz2<-runif(insp$myn)<insp$mypi
    mys2<-tbin + runif(insp$myn)*(tbin2-tbin) #switch in second period ##random switching
    mys2<-round(mys2, mypre) ##rounding
    
    mylarges2<-runif(insp$myn)>insp$lsprob
    mys2[mylarges2==TRUE]<-insp$lstime #large switch time
    
    myt2<-rexp(insp$myn, rate=((insp$mypar$theta *  (1-myz2+myz2*insp$mypar$tau) * insp$mypar$base))) #failure time, arm 2
    myt2<-round(myt2, mypre)
    myd2<-myt2<mymaxt2;    myt2x<-pmin(myt2,mymaxt2) ##censoring

    mydatad<-as.integer(c(myd1,myd2))
    mydatat<-c(myt1x, myt2x); mydatas<-c(mys1, mys2);  mydatar<-rep(1:2,each=insp$myn); mydatax<-rep(1,insp$myn*2)
    mydatag<-c((1+myz1)*(myt1x>=mys1), rep(0,insp$myn)) #0 = unknown, 1=don't switch, 2=switch
   
#reorder everything by time
    myord<-order(mydatat)
    mydatat<-mydatat[myord];  mydatas<-mydatas[myord]
    mydatax<-mydatax[myord];  mydatad<-mydatad[myord]
    mydatag<-mydatag[myord];  mydatar<-mydatar[myord]
    mynT<-2*insp$myn

    mydataz<-c(myz1, myz2)[myord]   
    
    myfitobj<-list(obs = list(mylik=as.double(-99), mydatat=mydatat, mydatad=mydatad, mydatax=mydatax, mydatag=mydatag, mydatar=mydatar, mynT=mynT, mym=mym, myp1=rep(0,insp$myn), mypitime=rep(0,insp$myn), myp3=rep(0,insp$myn), mypl=0, mydatas=mydatas), unobs=list(mydataz))

    myfitobj


  }



##used to generate stratified samples
fn.simdata.dbg<-function(insp)
  {
## debug
    ##    myn<-1000;  mypi<-0.25; mypar<-list(tau=0.5, theta=0.7, base=1)    
##    insp<-list(myn=1000, mypi=0.25, mypar=list(tau=0.5, theta=0.7, base=1), mycen=10)

    ## possible to change/add  later on
    tbin<-1; tbin2<-2;   mym<-1 #tbin2 = tbin then all switch at same tbin time
    
    ##censoring
#    mymaxt1<-insp$mycen*runif(insp$myn);  mymaxt2<-insp$mycen*runif(insp$myn)
    mymaxt1<-insp$mycen; mymaxt2<-insp$mycen

    ##switcher?
    myz1<-runif(insp$myn)<insp$mypi
    mys1<-tbin + runif(insp$myn)*(tbin2-tbin) #switch in second period ##random switching
    
    ##arm 1
    myt1as<-rexp(insp$myn, rate=( ( (1-myz1) + myz1*(insp$mypar$theta*insp$mypar$tau)) * insp$mypar$base))
    myt1ps<-rexp(insp$myn, rate=( ( (1-myz1) + myz1*insp$mypar$tau) * insp$mypar$base));  myt1ps<-pmin(myt1ps, mys1)
    myt1<-myt1ps;     myt1[myt1ps==mys1]<-(myt1ps+myt1as)[myt1ps==mys1]
    myd1<-myt1<mymaxt1 ;  myt1x<-pmin(myt1,mymaxt1)    ##censoring

    ##arm2 
##    myz2<-myz1;  mys2<-mys1     #switchers identical number to arm 1
    myz2<-runif(insp$myn)<insp$mypi;  mys2<-tbin + runif(insp$myn)*(tbin2-tbin) #switch in second period ##random switching   
    myt2<-rexp(insp$myn, rate=((insp$mypar$theta *  (1-myz2+myz2*insp$mypar$tau) * insp$mypar$base))) #failure time, arm 2 
    myd2<-myt2<mymaxt2;    myt2x<-pmin(myt2,mymaxt2) ##censoring

    mydatad<-as.integer(c(myd1,myd2))
    mydatat<-c(myt1x, myt2x); mydatas<-c(mys1, mys2);  mydatar<-rep(1:2,each=insp$myn); mydatax<-rep(1,insp$myn*2)
    mydatag<-c((1+myz1)*(myt1x>=mys1), rep(0,insp$myn)) #0 = unknown, 1=don't switch, 2=switch
   
#reorder everything by time
    myord<-order(mydatat)
    mydatat<-mydatat[myord];  mydatas<-mydatas[myord]
    mydatax<-mydatax[myord];  mydatad<-mydatad[myord]
    mydatag<-mydatag[myord];  mydatar<-mydatar[myord]
    mynT<-2*insp$myn

    mydataz<-c(myz1, myz2)[myord]   
    
    myfitobj<-list(obs = list(mylik=as.double(-99), mydatat=mydatat, mydatad=mydatad, mydatax=mydatax, mydatag=mydatag, mydatar=mydatar, mynT=mynT, mym=mym, myp1=rep(0,insp$myn), mypitime=rep(0,insp$myn), myp3=rep(0,insp$myn), mypl=0, mydatas=mydatas), unobs=list(mydataz))

    myfitobj


  }


##used to generate not stratified samples
## simulate data, IN: simulation pars (list), OUT: data object
fn.simdata<-function(insp)
  {
## debug
    ##    myn<-1000;  mypi<-0.25; mypar<-list(tau=0.5, theta=0.7, base=1)    
##    insp<-list(myn=1000, mypi=0.25, mypar=list(tau=0.5, theta=0.7, base=1), mycen=10)

    ## possible to add later on
    tbin<-1; tbin2<-2;   mym<-1 #tbin2 = tbin then all switch at same tbin time
    
    ##censoring
#    mymaxt1<-insp$mycen*runif(insp$myn);  mymaxt2<-insp$mycen*runif(insp$myn)
    mymaxt1<-insp$mycen; mymaxt2<-insp$mycen

    ##switcher?
    myz1<-runif(insp$myn)<insp$mypi
    mys1<-tbin + runif(insp$myn)*(tbin2-tbin) #switch in second period ##random switching
    
    ##arm 1
    myt1as<-rexp(insp$myn, rate=( ( (1-myz1) + myz1*(insp$mypar$theta*insp$mypar$tau)) * insp$mypar$base))
    myt1ps<-rexp(insp$myn, rate=( ( (1-myz1) + myz1*insp$mypar$tau) * insp$mypar$base));  myt1ps<-pmin(myt1ps, mys1)
    myt1<-myt1ps;     myt1[myt1ps==mys1]<-(myt1ps+myt1as)[myt1ps==mys1]
    myd1<-myt1<mymaxt1 ;  myt1x<-pmin(myt1,mymaxt1)    ##censoring

    ##arm2 

##    myz2<-myz1;  mys2<-mys1     #switchers identical number to arm 1
    myz2<-runif(insp$myn)<insp$mypi;  mys2<-tbin + runif(insp$myn)*(tbin2-tbin) #switch in second period ##random switching   
    myt2<-rexp(insp$myn, rate=((insp$mypar$theta *  (1-myz2+myz2*insp$mypar$tau) * insp$mypar$base))) #failure time, arm 2 
    myd2<-myt2<mymaxt2;    myt2x<-pmin(myt2,mymaxt2) ##censoring

    mydatad<-as.integer(c(myd1,myd2))
    mydatat<-c(myt1x, myt2x); mydatas<-c(mys1, mys2);  mydatar<-rep(1:2,each=insp$myn); mydatax<-rep(1,insp$myn*2)
    mydatag<-c((1+myz1)*(myt1x>=mys1), rep(0,insp$myn)) #0 = unknown, 1=don't switch, 2=switch
  
#reorder everything by time
    myord<-order(mydatat)
    mydatat<-mydatat[myord];  mydatas<-mydatas[myord]
    mydatax<-mydatax[myord];  mydatad<-mydatad[myord]
    mydatag<-mydatag[myord];  mydatar<-mydatar[myord]
    mynT<-2*insp$myn

    mydataz<-c(myz1, myz2)[myord]   
    
    myfitobj<-list(obs = list(mylik=as.double(-99), mydatat=mydatat, mydatad=mydatad, mydatax=mydatax, mydatag=mydatag, mydatar=mydatar, mynT=mynT, mym=mym, myp1=rep(0,insp$myn), mypitime=rep(0,insp$myn), myp3=rep(0,insp$myn), mypl=0, mydatas=mydatas), unobs=list(mydataz))

    myfitobj


  }




################
################

#####################################
## partial likelihood simulation, main routine. IN: seed, simulation setup, OUT: estimates
fn.plsim<-function(inseed, insp)
  {

##   inseed<-15;     insp<-list(myn=1000, mypi=0.25, mypar=list(tau=0.5, theta=0.7, base=.1), mycen=3) ##debug
##    insp<-insp1[[3]]
    print(paste("RSEED",inseed))
    set.seed(inseed)
    
   sim.dta.all<-fn.simdata(insp)
   sim.dta<-sim.dta.all[[1]]
   sim.z<-sim.dta.all[[2]] 

    ##starting values, simple estimates
    
    ##itt effect
    init.theta<-(sum(sim.dta$mydatad[sim.dta$mydatar==2]) / sum(sim.dta$mydatat[sim.dta$mydatar==2])) / (sum(sim.dta$mydatad[sim.dta$mydatar==1]) / sum(sim.dta$mydatat[sim.dta$mydatar==1]))

    ## effect after switching in switchers vs everyone in treatment
    init.switch<-(sum(sim.dta$mydatad[sim.dta$mydatag==2 & sim.dta$mydatat>sim.dta$mydatas]) / sum(sim.dta$mydatat[sim.dta$mydatag==2 & sim.dta$mydatat>sim.dta$mydatas])) /  (sum(sim.dta$mydatad[sim.dta$mydatar==2 & sim.dta$mydatat>sim.dta$mydatas]) / sum(sim.dta$mydatat[sim.dta$mydatar==2 & sim.dta$mydatat>sim.dta$mydatas]))

    if(init.switch==0) init.switch<-0.1

    ##temp.pl1<-fn.PL(c(log(init.switch),log(init.theta)), sim.dta )
##    if(init.switch>1) {
 ##     mypstart<-log(c(init.switch^4,init.theta))
  ##  } else {
   ##   mypstart<-log(c(1,init.theta))
    ##}
    
    mypstart<-log(c(init.switch,init.theta))
    
    temp.opt<-fn.PL.fit(mypstart, sim.dta )
##    temp.opt<-fn.PL.fit(c(0,0), sim.dta )
    temp.pl<- fn.PL.ALL(temp.opt[[1]], sim.dta)

    
    ## if using quasi-Newon (no diff really)
##    temp.opt<-optim(log(c(init.switch^4, init.theta)), fn=fn.PL,  method="L-BFGS", ind=sim.dta, lower=log(c(0.0001, 0.0001)), upper=log(c(100,100)))
##        temp.opt<-optim(temp.opt[[1]], fn=fn.PL,  method="L-BFGS", ind=sim.dta, lower=log(c(0.0001, 0.0001)), upper=log(c(100,100)))
  ##  temp.pl<- fn.PL.ALL(temp.opt[[1]], sim.dta)
##    temp.pl<- fn.PL.ALL(c(0,0), sim.dta)
##    c(temp.opt[[1]], myopt$par, temp.pl[[11]][1], sum(sim.z[[1]][sim.dta$mydatar==1]), sum(sim.z[[1]][sim.dta$mydatar==2]))

    c(temp.opt[[1]], temp.pl[[11]][1], sum(sim.z[[1]][sim.dta$mydatar==1]), sum(sim.z[[1]][sim.dta$mydatar==2]))
    
  }


## partial likelihood simulation, main routine with bootstrap. IN: seed, simulation setup, OUT: estimates and se's
fn.plsim.boot<-function(inseed, insp, inR=99)
  {

##   inseed<-15;     insp<-list(myn=1000, mypi=0.25, mypar=list(tau=0.5, theta=0.7, base=.1), mycen=3) ##debug
##    insp<-insp1[[1]]
    print(paste("RSEED",inseed))
    set.seed(inseed)
    
   sim.dta.all<-fn.simdata(insp)
   sim.dta<-sim.dta.all[[1]]
   sim.z<-sim.dta.all[[2]] 

    ##starting values, simple estimates
    
    ##itt effect
    init.theta<-(sum(sim.dta$mydatad[sim.dta$mydatar==2]) / sum(sim.dta$mydatat[sim.dta$mydatar==2])) / (sum(sim.dta$mydatad[sim.dta$mydatar==1]) / sum(sim.dta$mydatat[sim.dta$mydatar==1]))

    ## effect after switching in switchers vs everyone in treatment
    init.switch<-(sum(sim.dta$mydatad[sim.dta$mydatag==2 & sim.dta$mydatat>sim.dta$mydatas]) / sum(sim.dta$mydatat[sim.dta$mydatag==2 & sim.dta$mydatat>sim.dta$mydatas])) /  (sum(sim.dta$mydatad[sim.dta$mydatar==2 & sim.dta$mydatat>sim.dta$mydatas]) / sum(sim.dta$mydatat[sim.dta$mydatar==2 & sim.dta$mydatat>sim.dta$mydatas]))

    if(init.switch==0) init.switch<-0.1

    ##if(init.switch>1) {
     ## mypstart<-log(c(init.switch^4,init.theta))
    ##} else {
     ## mypstart<-log(c(1,init.theta))
   ## }
    ##temp.pl1<-fn.PL(c(log(init.switch),log(init.theta)), sim.dta )

    mypstart<-log(c(init.switch,init.theta))
    
    temp.opt<-fn.PL.fit(mypstart, sim.dta )
    temp.pl<- fn.PL.ALL(temp.opt[[1]], sim.dta)


    fn.fitmeboot<-function(myrids, inidx, sim.dta, myfit)
      {

##        inidx<-seq(1,2000); myfit<-temp.opt[[1]] ##debug; ##        inidx<-c(ceiling(runif(1000)*1000), ceiling(1000+runif(1000)*1000))
       
        tidx<-myrids[inidx]       
        
        ## create bootstrap data
        mydta<-sim.dta

        #time
        mydta$mydatat<-sim.dta$mydatat[tidx];
        myord<-order(mydta$mydatat)
        mydta$mydatat<-mydta$mydatat[myord]

        #rest
        mydta$mydatad<-sim.dta$mydatad[tidx[myord]]
        mydta$mydatax<-sim.dta$mydatax[tidx[myord]]; mydta$mydatar<-sim.dta$mydatar[tidx[myord]]
        mydta$mydatas<-sim.dta$mydatas[tidx[myord]]; mydta$mydatag<-sim.dta$mydatag[tidx[myord]]

        #fit
  ##      temp.pl<- fn.PL.ALL(myfit, mydta)
    ##    temp.pl
        temp.opt.boot<-fn.PL.fit(myfit, mydta) ##internal alg
        temp.pl.boot<- fn.PL.ALL(temp.opt.boot[[1]], mydta)

        c(temp.opt.boot[[1]], temp.pl.boot[[11]][1])

      }
    
    myrids<-c(seq(1,sim.dta$myn)[sim.dta$mydatar==1], seq(1,sim.dta$myn)[sim.dta$mydatar==2])
    mystrata<-rep(1:2,each=sim.dta$myn/2)

    t.boot<-boot(data=myrids, statistic=fn.fitmeboot, strata=mystrata, sim.dta=sim.dta, myfit=temp.opt[[1]], R=inR)
    outvar<-apply(t.boot$t, 2, var)
    outpc<-sapply(1:3, function(idx) boot.ci(t.boot, index=idx, type="perc")$percent[4:5])

    ## fits, real number, bootstrap variance and percentile ci's
    c(temp.opt[[1]], temp.pl[[11]][1], sum(sim.z[[1]][sim.dta$mydatar==1]), sum(sim.z[[1]][sim.dta$mydatar==2]), outvar, outpc)
    
  }



##partial likelihood fit, IN: parameters, data; OUT: fit 
fn.PL<-function(inpar, ind)
      {

        temp<-.C("rint_pl", as.double(inpar), as.double(ind$mylik), as.double(ind$mydatat), as.integer(ind$mydatad), as.double(ind$mydatax), as.integer(ind$mydatag), as.integer(ind$mydatar), as.integer(ind$mynT), as.integer(ind$mym), as.double(ind$mydatas), as.double(ind$myp1), as.double(ind$mypitime), as.double(ind$mypl))
        print(paste("PL",round(temp[[13]],2), "pars:", exp(inpar[1]), exp(inpar[2])))
        temp[[13]]
      }
#partial likelihood, all the object
fn.PL.ALL<-function(inpar, ind)
      {

        temp<-.C("rint_pl", as.double(inpar), as.double(ind$mylik), as.double(ind$mydatat), as.integer(ind$mydatad), as.double(ind$mydatax), as.integer(ind$mydatag), as.integer(ind$mydatar), as.integer(ind$mynT), as.integer(ind$mym), as.double(ind$mydatas), as.double(ind$myp1), as.double(ind$mypitime), as.double(ind$mypl))
        temp
      }


##partial likelihood fit, IN: parameters, data; OUT: fit 
fn.PL.fit<-function(inpar, ind)
      {

        temp<-.C("testerpl", as.double(inpar), as.double(ind$mylik), as.double(ind$mydatat), as.integer(ind$mydatad), as.double(ind$mydatax), as.integer(ind$mydatag), as.integer(ind$mydatar), as.integer(ind$mynT), as.integer(ind$mym), as.double(ind$mydatas), as.double(ind$myp1), as.double(ind$mypitime), as.double(ind$mypl))
     ##   print(paste("PL",round(temp[[13]],2), "pars:", exp(inpar[1]), exp(inpar[2])))
        temp
      }

## FULL LIKELIHOOD Sim
fn.fullsim<-function(inseed, insp)
  {

##   inseed<-15;     insp<-list(myn=1000, mypi=0.25, mypar=list(tau=0.5, theta=0.7, base=.1), mycen=3) ##debug
##    insp<-insp1[[3]]; rseed<-1002
 #   insp$mypar$base<-1
#    insp$mycen<-3
    print(paste("RSEED",inseed))
    set.seed(inseed)
    
   sim.dta.all<-fn.simdata(insp)
   sim.dta<-sim.dta.all[[1]]
   sim.z<-sim.dta.all[[2]] 

   sim.dta$mydatasidx<-sapply(1:sim.dta$myn, function(idx) pmin(sim.dta$myn-1,sum(sim.dta$mydatat<sim.dta$mydatas[idx])))
##   sim.dta$eventlist<-seq(1,length(sim.dta$mydatad))[sim.dta$mydatad==1]
 ##  sim.dta$nevents<-length(sim.dta$eventlist)
    
    
    ##starting values, simple estimates
    
    ##itt effect
    init.theta<-(sum(sim.dta$mydatad[sim.dta$mydatar==2]) / sum(sim.dta$mydatat[sim.dta$mydatar==2])) / (sum(sim.dta$mydatad[sim.dta$mydatar==1]) / sum(sim.dta$mydatat[sim.dta$mydatar==1]))

    ## effect after switching in switchers vs everyone in treatment
    init.switch<-(sum(sim.dta$mydatad[sim.dta$mydatag==2 & sim.dta$mydatat>sim.dta$mydatas]) / sum(sim.dta$mydatat[sim.dta$mydatag==2 & sim.dta$mydatat>sim.dta$mydatas])) /  (sum(sim.dta$mydatad[sim.dta$mydatar==2 & sim.dta$mydatat>sim.dta$mydatas]) / sum(sim.dta$mydatat[sim.dta$mydatar==2 & sim.dta$mydatat>sim.dta$mydatas]))

    if(init.switch==0) init.switch<-0.1

    ##temp.pl1<-fn.PL(c(log(init.switch),log(init.theta)), sim.dta )

    ##starting vals
##    if(init.switch>1) {
 ##     mypstart<-log(c(init.switch^4,init.theta))
  ##  } else {
   ##   mypstart<-log(c(1,init.theta))
   ## }
#        mypstart<-log(c(0.2,0.7))
#        mypstart<-log(c(0.01,0.7))
#    mypstart<-log(mypstart/100 / (1-mypstart/100))
#   temp<-fn.FULL(mypstart, sim.dta)
    mypstart<-log(c(init.switch,init.theta)   )
    
    myopt<-optim(c(0,0), fn=fn.FULL,  method="BFGS", ind=sim.dta)
##        myopt<-optim(mypstart, fn=fn.FULL,  method="Nelder-Mead", ind=sim.dta)
##    myopt<-optim(mypstart, fn=fn.FULL,  method="BFGS", ind=sim.dta)
##   temp<-fn.FULL.all(mypstart, sim.dta)
 ##     temp<-fn.FULL.all(c(1,3), sim.dta)
    
   temp<-fn.FULL.all(myopt$par, sim.dta)    

    c(myopt$par, temp[[11]][1], sum(sim.z[[1]][sim.dta$mydatar==1]), sum(sim.z[[1]][sim.dta$mydatar==2]))
    
  }






#######################




