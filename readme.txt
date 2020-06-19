This program requires rstan, foreach and doparallel
The program is launched in BHM.R, which pre-processes the data
which then calls BHM.stan, which runs the Bayesian Hierarchical Mixture (BHM).

**Note that for the first run, the STAN program needs to be compiled. Make sure you set initiate_stan (in terminal type: init_stan=T) before running your R script. Then turn off afterwards (in terminal type: init_stan=F). This step does not need to be redone, unless the STAN code is changed.
**also adjust the number of cores (options(cores=30)) in the R script to be less than the total number on your computer. 

The publicly data can be downloaded from the Living Planet Index database www.livingplanetindex.org/. The population IDs used in this analyses can be found in the file  “info.rds”, which can be matched to the LPI database. To run the program, raw growth rates need to be calculated for each population, and placed into a list, where the names of each element in the list need to correspond to the IDs from the LPI dataset (see “growth_rate_example.rds”, for the structure of the list). The code to generate the growth rate file is:

		growth_rate=by(tot,tot$id,function(x)
		{
				if(length(which(x[,cN]==0))>0) #if any measurement is zero, offset the entire time series by 1% of the mean
				{
					x[,cN]=x[,cN]+0.01*mean(x[,cN])
				}
			}
			x=x[order(x$Year),] #make sure data are in temporal order
			dur=x$Year[-1]-x$Year[-nrow(x)] #when years are missing, take the mean growth between those intervals
			gr<-log(x[-1,cN]/x[-nrow(x),cN])/dur
		return(gr)})

where tot is the LPI dataset, tot$id is the population ids, cN is the column number of the abundance estimates, and Year is the sample year.


After STAN runs, the output is placed in a summary file called “dat.rds”. Users should check the Rhat values. if Rhat deviates substantially from one, convergence is suspect. Please re-run those individual analyses.


For questions, please contact: Brian Leung (brian.leung2@mcgill.ca)
