coxsc Copyright (C) 2019 Adam Brentnall

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


Summary of project

R code to fit a binary and proportional hazards model for 2-arm trials with selective crossover

Still working towards a proper R package with everything (but any help welcome!)

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To compile the code (for proportional hazards model) run the commands in installation.txt

Then to use it go to the installed package in the subdirectory you created (built).

Or use the precompiled code from my machine in <built-mymachine>

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
R code to help in the demo folder. This has the following. 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A. Binary model and example with plots etc

* 180413-binary-xover.r - R functions and analysis code. Does not need any C functions etc.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
B. Proportional hazard model simulations in paper 

* pifuncs.R  - R functions
* Rsim_rmpi.r  - Call this using open MPI (orterun, mpirun, mpiexec). Will need it installed on computer, and Rmpi. It will take a long time so there are options for which simulations you would like to run. e.g. Only partial likelihood (quickest!)

Example of how to run:

> orterun -n 1 -H localhost,localhost,localhost,localhost,localhost,localhost,localhost,localhost,localhost,localhost  R --slave -f Rsim_rmpi.r

* Rsim_analysis.r  - Code to do analysis of the simulation results, same format as paper

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
B. Some further code developments for PH model, to fit stratified baseline hazard function, profile likelihood confidence intervals. Still in progress, but appears to work OK

* fitfns.r  - R functions
* strat-pl-eg.r - Example to fit a model
