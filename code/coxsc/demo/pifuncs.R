#############################################
## Some code used in the paper
## Adam Brentnall
## 12th April 2018
#############################################


## Used by profile likelihood. Treatment effect fixed, pass in ind as $trteff
## Stratified baseline
fn.FULL.trtfixed<-function(inpar, ind)
      {

        fitpar<-c(inpar[1], ind$trteff) ##no covariates yet
            
        temp<-.C("testerfull_strata", as.double(fitpar), as.double(ind$mylik), as.double(ind$mydatat), as.integer(ind$mydatad), as.double(ind$mydatax), as.integer(ind$mydatag), as.integer(ind$mydatar), as.integer(ind$mynT), as.integer(ind$mym), as.double(ind$mydatas), as.double(ind$myp1), as.double(ind$mypitime), as.double(ind$mypl), as.integer(ind$mydatasidx), as.double(ind$myp3), as.integer(ind$mystrata), as.integer(ind$mynstrata))
        print(paste("FULL LIK",round(temp[[2]],2), "pars:",exp(fitpar[1]), exp(fitpar[2])))
       temp[[2]]
      }

## Used by profile likelihood. Treatment effect fixed, pass in ind as $trteff
## Not stratified baseline
fn.FULLxs.trtfixed<-function(inpar, ind)
      {

        fitpar<-c(inpar[1], ind$trteff) ##no covariates yet
            
        temp<-.C("testerfull", as.double(fitpar), as.double(ind$mylik), as.double(ind$mydatat), as.integer(ind$mydatad), as.double(ind$mydatax), as.integer(ind$mydatag), as.integer(ind$mydatar), as.integer(ind$mynT), as.integer(ind$mym), as.double(ind$mydatas), as.double(ind$myp1), as.double(ind$mypitime), as.double(ind$mypl), as.integer(ind$mydatasidx), as.double(ind$myp3), as.integer(ind$mystrata), as.integer(ind$mynstrata))
        print(paste("FULL LIK",round(temp[[2]],2), "pars:",exp(fitpar[1]), exp(fitpar[2])))
       temp[[2]]
      }


## profile likelihood confidence interval
## mlpar is MLE all pars, mllik is corresponding log likelihood, inalp is for confint
## likelihood fit, IN: parameters, data; OUT: fit 
fn.FULL.proflik<-function(mlpar, mllik, ind, inalp=0.05, coxstrata=TRUE)
      {

        ##debug --
#        inalp<-0.05; ind<-sim.dta       
        ##mllik<-myopt.st$value
 #       mllik<-myopt.st.obj[[2]]
  #      mlpar<-myopt.st.obj[[1]]
 
#######################
        ind$trteff<-mlpar[2]
        mychi<-qchisq(1-inalp, 1) 

        fn.rootfind.1<-function(intrt, inother, mllik, mychi, ind, coxstrata=TRUE)
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

       mylower<-log(0.1) #lowest possible -- to investigate
        myupper<-log(10) #lowest possible -- to investigate changing 
      
  ##        temp<- fn.rootfind.1( mlpar[2], mlpar[1], mllik, mychi, ind)
       ci.lower<-uniroot(fn.rootfind.1, interval=c(mylower, mlpar[2]) , inother=mlpar[1], mllik=mllik, mychi=mychi, ind=ind)
       ci.upper<-uniroot(fn.rootfind.1, interval=c(mlpar[2], myupper) , inother=mlpar[1], mllik=mllik, mychi=mychi, ind=ind)

       myout<-(c(ci.lower[[1]], mlpar[2],ci.upper[[1]]))
        myout
      }

## likelihood fit, IN: parameters, data; OUT: fit 
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
##        print(paste("FULL LIK",round(temp[[2]],2), "pars:",exp(inpar)))
        temp
      }

##write simulated data to output for input to main prog
fn.writecsvinout<-function(ind)
  {
ind<-sim.dta ##debug
     myout<-data.frame(id=seq(1,ind$mynT), 
     t=ind$mydatat,
     d=ind$mydatad,
     g=ind$mydatag,
     r=ind$mydatar,
     sx=ind$mydatas,
     st=rep(seq(1,ind$mynstrata), diff(ind$mystrata)))

   write.csv(myout, file="demoin.csv", row.names=FALSE)

  }

## partial likelihood simulation, main routine. IN: seed, simulation setup, OUT: estimates
fn.fullsim.dbg<-function(inseed, insp)
  {
   print(paste("RSEED",inseed))
    set.seed(inseed)
## insp<-insp1[[8]] ##debug
   
   sim.dta.all1<-fn.simdata.dbg(insp)
   ##add strat1a
   insp$mypar$base<-insp$mypar$base*2
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
   

    
    ##starting values, simple estimates
    
    ##itt effect
    init.theta<-(sum(sim.dta$mydatad[sim.dta$mydatar==2]) / sum(sim.dta$mydatat[sim.dta$mydatar==2])) / (sum(sim.dta$mydatad[sim.dta$mydatar==1]) / sum(sim.dta$mydatat[sim.dta$mydatar==1]))

    ## effect after switching in switchers vs everyone in treatment
    init.switch<-(sum(sim.dta$mydatad[sim.dta$mydatag==2 & sim.dta$mydatat>sim.dta$mydatas]) / sum(sim.dta$mydatat[sim.dta$mydatag==2 & sim.dta$mydatat>sim.dta$mydatas])) /  (sum(sim.dta$mydatad[sim.dta$mydatar==2 & sim.dta$mydatat>sim.dta$mydatas]) / sum(sim.dta$mydatat[sim.dta$mydatar==2 & sim.dta$mydatat>sim.dta$mydatas]))

    if(init.switch==0) init.switch<-0.1

    mypstart<-log(c(init.switch,init.theta)   )


#    temp<-fn.FULL.all(myopt$par, sim.dta)    

##stratified#####################
   sim.dta<-sim.dta.all1[[1]]
   mydatasidx1<-sapply(1:sim.dta$myn, function(idx) pmin(sim.dta$myn-1,sum(sim.dta$mydatat<sim.dta$mydatas[idx])))
   sim.dta$mydatasidx<-mydatasidx1
  ##  myopt<-optim(mypstart, fn=fn.FULL,  method="BFGS", ind=sim.dta) ##fit

   sim.dta<-sim.dta.all2[[1]]
   mydatasidx2<-sapply(1:sim.dta$myn, function(idx) pmin(sim.dta$myn-1,sum(sim.dta$mydatat<sim.dta$mydatas[idx])))
  sim.dta$mydatasidx<-mydatasidx2
##   myopt<-optim(mypstart, fn=fn.FULL,  method="BFGS", ind=sim.dta) ##fit
   
   sim.dta<-sim.dta.all[[1]]
##   sim.dta$mydatasidx<-c(mydatasidx1, mydatasidx2)
   sim.dta$mystrata<-c(0,sim.dta$myn/2, sim.dta$myn)
   sim.dta$mynstrata<-2

    myopt.st.pl<- fn.PL.fit.dbg(mypstart, sim.dta)
   mypstart.pl<-myopt.st.pl[[1]]
    myopt.st<-optim(mypstart.pl, fn=fn.FULL.dbg,  method="BFGS", ind=sim.dta) ##use PL as start
##    myopt.st<-optim(mypstart, fn=fn.FULL.dbg,  method="BFGS", ind=sim.dta) ##use other as start  
    myopt.st.pl<- fn.PL.fit.dbg(mypstart, sim.dta)
    myopt.st.obj<- fn.FULL.all.dbg(myopt.st$par, sim.dta)    

##########bootstrap PL variance etc
    fn.fitmeboot<-function(myrids, inidx, sim.dta, myfit)
      {

        tidx<-myrids[inidx]       
        
        ## create bootstrap data
        mydta<-sim.dta

        #time
        mydta$mydatat<-sim.dta$mydatat[tidx];
        myord1<-order(mydta$mydatat[1:(sim.dta$mystrata[2])])
        myord2<-order(mydta$mydatat[(sim.dta$mystrata[2]+1) : (sim.dta$mystrata[3])])
        myord<-c(myord1, sim.dta$mystrat[2]+myord2)
        
#        myord<-order(mydta$mydatat)
        
        mydta$mydatat<-mydta$mydatat[myord]

        #rest
        mydta$mydatad<-sim.dta$mydatad[tidx[myord]]
        mydta$mydatax<-sim.dta$mydatax[tidx[myord]]; mydta$mydatar<-sim.dta$mydatar[tidx[myord]]
        mydta$mydatas<-sim.dta$mydatas[tidx[myord]]; mydta$mydatag<-sim.dta$mydatag[tidx[myord]]

        #fit
  ##      temp.pl<- fn.PL.ALL(myfit, mydta)
    ##    temp.pl
##        temp.opt.boot<-fn.PL.fit(myfit, mydta) ##internal alg
##        temp.pl.boot<- fn.PL.ALL(temp.opt.boot[[1]], mydta)
        temp.pl.boot<- fn.PL.fit.dbg(myfit, mydta)

       temp.pl.boot[[1]]

      }

##    myrids<-c(seq(1,sim.dta$mynT)[sim.dta$mydatar==1], seq(1,sim.dta$mynT)[sim.dta$mydatar==2])
    myrids<-c(seq(1,sim.dta$mynT))
    mystrata<-2*(sim.dta$mydatar-1) + rep(c(1,2), each=2000)

    temp<-fn.fitmeboot(myrids, seq(1,4000), sim.dta, myopt.st.pl[[1]])

   inR<-99    
    t.boot<-boot(data=myrids, statistic=fn.fitmeboot, strata=mystrata, sim.dta=sim.dta, myfit=myopt.st.pl[[1]], R=inR)
    outvar<-apply(t.boot$t, 2, var)
    outpc<-sapply(1:2, function(idx) boot.ci(t.boot, index=idx, type="perc")$percent[4:5])

    outproflikci<-fn.FULL.proflik(myopt.st.obj[[1]], myopt.st.obj[[2]], sim.dta, inalp=0.05)

##
     ## not stratified
    sim.dta$mydatasidx<-sapply(1:sim.dta$myn, function(idx) pmin(sim.dta$myn-1,sum(sim.dta$mydatat<sim.dta$mydatas[idx])))
   myopt<-optim(mypstart, fn=fn.FULL,  method="BFGS", ind=sim.dta)
   myopt.obj<- fn.FULL.all(myopt$par, sim.dta)    
   myopt.pl<-fn.PL.fit(mypstart, sim.dta )
   

    c(exp(myopt$par), myopt.st.obj[[11]][1], exp(myopt.st$par), myopt.st.obj[[11]][1],
      exp(myopt.pl[[1]]), myopt.pl[[11]][1], myopt.pl[[11]][1001],
      exp(myopt.st.pl[[1]]), (myopt.st.pl[[11]][1]), myopt.st.pl[[11]][1001],
      sum(sim.z[[1]][sim.dta$mydatar==1]), sum(sim.z[[1]][sim.dta$mydatar==2]))
   
    
  }





                                        #partial likelihood, all the object
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



## partial likelihood simulation, main routine. IN: seed, simulation setup, OUT: estimates
fn.plsim.dbg<-function(inseed, insp)
  {

##   inseed<-15;     insp<-list(myn=1000, mypi=0.25, mypar=list(tau=0.5, theta=0.7, base=.1), mycen=3) ##debug
##    insp<-insp1[[8]]
   print(paste("RSEED",inseed))
    set.seed(inseed)
    
   sim.dta.all1<-fn.simdata.dbg(insp)
   ##add strat1a
   insp$mypar$base<-insp$mypar$base*2
   sim.dta.all2<-fn.simdata.dbg(insp)



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
                       mydatas=c(sim.dta.all1$obs$mydatas, sim.dta.all2$obs$mydatas)),
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
                       mydatas=c(sim.dta.all1$obs$mydatas, sim.dta.all2$obs$mydatas))[myord],
                     unobs=list(c(unlist(sim.dta.all1$unobs), unlist(sim.dta.all2$unobs))[myord] ))

   
   sim.dta<-sim.dta.all.sort[[1]]
   sim.z<-sim.dta.all.sort[[2]] 

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
##     temp.opt1<-fn.PL.fit(mypstart, sim.dta.all1[[1]] )
##     temp.opt2<-fn.PL.fit(mypstart, sim.dta.all2[[1]] )
##    temp.opt<-fn.PL.fit(c(0,0), sim.dta )

   sim.dta<-sim.dta.all[[1]]
   sim.dta$mystrata<-c(0,sim.dta$mynT/2, sim.dta$mynT)
   sim.dta$mynstrata<-2


##    temp.pl<- fn.PL.ALL.dbg(mypstart, sim.dta)
     temp.pl<- fn.PL.fit.dbg(mypstart, sim.dta)
#   temp.p2<- fn.PL.ALL.dbg(temp.opt1[[1]], sim.dta)
#   temp.o2<-fn.PL.ALL(temp.opt1[[1]], sim.dta.all1[[1]] )
 #  temp.pl2<- fn.PL.ALL.dbg(temp.opt[[1]], sim.dta)
    
    ## if using quasi-Newon (no diff really)
##    temp.opt<-optim(log(c(init.switch^4, init.theta)), fn=fn.PL,  method="L-BFGS", ind=sim.dta, lower=log(c(0.0001, 0.0001)), upper=log(c(100,100)))
##        temp.opt<-optim(temp.opt[[1]], fn=fn.PL,  method="L-BFGS", ind=sim.dta, lower=log(c(0.0001, 0.0001)), upper=log(c(100,100)))
  ##  temp.pl<- fn.PL.ALL(temp.opt[[1]], sim.dta)
##    temp.pl<- fn.PL.ALL(c(0,0), sim.dta)
##    c(temp.opt[[1]], myopt$par, temp.pl[[11]][1], sum(sim.z[[1]][sim.dta$mydatar==1]), sum(sim.z[[1]][sim.dta$mydatar==2]))

#   temp.opt[[11]][1]
#   temp.opt1[[11]][1]
#   temp.pl[[11]][1]
#   temp.opt2[[11]][1]
#   temp.pl[[11]][1001]   
   
    c(temp.opt[[1]], temp.opt[[11]][1], temp.pl[[1]], temp.pl[[11]][1], temp.pl[[11]][1001], sum(sim.z[[1]][sim.dta$mydatar==1]), sum(sim.z[[1]][sim.dta$mydatar==2]))
    
  }

fn.simdata.dbg<-function(insp)
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



fn.simdata.dbg<-function(insp)
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





#############################################################
## BELOW PRE PACKAGE DEVELOPMENT 22/09/2013
#############################################################
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


#####################################
## FULL LIKELIHOOD TESTING
#####################################

##partial likelihood fit, IN: parameters, data; OUT: fit 
fn.FULL<-function(inpar, ind)
      {
                   
            temp<-.C("testerfull", as.double(inpar), as.double(ind$mylik), as.double(ind$mydatat), as.integer(ind$mydatad), as.double(ind$mydatax), as.integer(ind$mydatag), as.integer(ind$mydatar), as.integer(ind$mynT), as.integer(ind$mym), as.double(ind$mydatas), as.double(ind$myp1), as.double(ind$mypitime), as.double(ind$mypl), as.integer(ind$mydatasidx), as.double(ind$myp3))
            print(paste("FULL LIK",round(temp[[2]],2), "pars:",exp(inpar[1]), exp(inpar[2])))
            temp[[2]]
      }

fn.FULL.all<-function(inpar, ind)
      {

#        thispar<-log(exp(inpar)/(1+exp(inpar)) * 100)
        temp<-.C("testerfull", as.double(inpar), as.double(ind$mylik), as.double(ind$mydatat), as.integer(ind$mydatad), as.double(ind$mydatax), as.integer(ind$mydatag), as.integer(ind$mydatar), as.integer(ind$mynT), as.integer(ind$mym), as.double(ind$mydatas), as.double(ind$myp1), as.double(ind$mypitime), as.double(ind$mypl), as.integer(ind$mydatasidx), as.double(ind$myp3))
##        print(paste("FULL LIK",round(temp[[2]],2), "pars:",exp(inpar)))
        temp
      }


## partial likelihood simulation, main routine. IN: seed, simulation setup, OUT: estimates
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






