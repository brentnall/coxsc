## Do simulations reported in the paper

library("Rmpi")

##use 10 cores
mpi.spawn.Rslaves(nslaves=10)

if (mpi.comm.size() < 2) {
      print("More slave processes are required.")
          mpi.quit()
    }

.Last <- function(){
      if (is.loaded("mpi_initialize")){
                if (mpi.comm.size(1) > 0){
                              print("Please use mpi.close.Rslaves() to close slaves.")
                                          mpi.close.Rslaves()
                            }
                        print("Please use mpi.quit() to quit R")
                        .Call("mpi_finalize")
              }
    }

##libraries
source("pifuncs.R"); dyn.load("../libs/coxsc.so"); library("boot")

doBIG<-FALSE ## do large sample sims?

doREP<-FALSE## do PL sims?

FL<-TRUE ##do full lik sims?

nreps<-1000

insp1<-fn.simscen(1000); insp2<-fn.simscen(100000)

##send everything to the slaves

mpi.bcast.Robj2slave(insp1); mpi.bcast.Robj2slave(nreps); mpi.bcast.Robj2slave(insp2)

mpi.bcast.cmd(source("pifuncs.R")); mpi.bcast.cmd(dyn.load("../libs/coxsc.so"));  mpi.bcast.cmd(library("boot"))

## now do the simulations

myout1<-list(); myout2<-list()

if(FL==TRUE) ## Full likelihood run
  {
   for(idx in 1:length(insp1))
     {
       print(paste("FL Starting run",idx,"/", length(insp1), "..."))

       mpi.bcast.Robj2slave(idx)

       myout1[[idx]]<-mpi.applyLB(1:nreps, function(idy) fn.fullsim(idx*nreps+idy, insp1[[idx]])) 

       save(file="simresFL_scen222.Rdata",myout1)       

     }
  }

    if(doREP){ ## Partial likelihood run
      for(idx in 1:length(insp1))
        {
          print(paste("PL Starting run",idx,"/", length(insp1), "..."))

          mpi.bcast.Robj2slave(idx)


          myout1[[idx]]<-mpi.applyLB(1:nreps, function(idy) fn.plsim(idx*nreps+idy, insp1[[idx]]))

          save(file="simres_scen22.Rdata",myout1)
          
          ##myout1[[idx]]<-mpi.applyLB(1:nreps, function(idy) fn.plsim.boot(idx*nreps+idy, insp1[[idx]],199)) ##bootstrap ci as well

          ##save(file="simresBOOT_scen22.Rdata",myout1)                    

        }
    }
    
    if(doBIG){ ## large sample
      myout2<-mpi.applyLB(1:length(insp2), function(idy) fn.plsim(idy, insp2[[idy]]))

      save(file="simres2_scen2.Rdata",myout2)
    }
  

##tidy up
mpi.close.Rslaves()

mpi.quit(save="no")
