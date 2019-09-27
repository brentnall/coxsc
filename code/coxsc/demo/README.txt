R code to demonstrate fitting the proportional hazards model for 2-arm trials with selective crossover
Adam Brentnall
13th April 2018

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A. To repeat simulations in paper with code used

* pifuncs.R  - R functions
* Rsim_rmpi.r  - Call this using open MPI (orterun, mpirun, mpiexec). Will need it installed on computer, and Rmpi. It will take a long time so there are options for which simulations you would like to run. e.g. Only partial likelihood (quickest!)
* runme.txt  - Example of how to call the simulations on 10 cores (command from bash)
* Rsim_analysis.r  - Code to do analysis of the simulation results, same format as paper

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
B. Some further code developments, to fit stratified baseline hazard function, profile likelihood confidence intervals. 

* fitfns.r  - R functions
* strat-pl-eg.r - Example to fit a model
