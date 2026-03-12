#####################################################################################################################
# Novel g-computation algorithms for time-varying actions with recurrent and semi-competing events
# D'Alessio et al.
# March 1, 2026
#####################################################################################################################

library(tidyverse)
library(resample)
library(nnet)
library(stats)

setwd("YOUR DIRECTORY HERE")

source("semicomp_gcomp_estimator_functions.R")
source("sim_functions.R")

#Set seed for reproducibility
set.seed(10312025)

B <- 500 #500 bootstraps
time<-3 #3 follow-up times after baseline
estimators<-c("standard","ice") #list of estimators #"naive","ice",
state_count<-3 #number of states (No hypertension, Hypertension, Death)
event<-c(1,2,3) #number corresponding to (intermediate event, terminal event)
plan_list<-c("c(0,0,0)","c(1,1,1)") #untreated, treated plan
prefixes <- c("TrNon","TrInt", "TrTer", "NtrNon","NtrInt", "NtrTer")  #TrNon = treated, no intermediate or terminal event; TrInt = treated, intermediate event; TrTer = treated, terminal; NtrNon = Not treated, no intermediate or terminal event; NtrInt = Not treated, intermediate event; NtrTer = Not treated, terminal

#Create empty data frame for the bootstrap results
ice.bsresults      <- matrix(nrow = B, ncol = time*state_count*2) # ICE g-computation estimator (we multiply by 2 because we calculate mean and se for each measure)
standard.bsresults <- matrix(nrow = B, ncol = time*state_count*2) # Standard g-computation estimator
  
#Name columns of bootstrap results (not used but helps when checking mid-steps)
result_matrix_names <- unlist(lapply(prefixes, function(prefix) {
    unlist(lapply(1:time, function(i) {
      paste0(prefix, "_T", i) #format of TH_T1, TH_T2, etc.
    }))
}))
  
dimnames(ice.bsresults) <- list(NULL, result_matrix_names)
dimnames(standard.bsresults) <- list(NULL, result_matrix_names)
  
#Create data frame for final results
ice_results<-create_results_df(time)
standard_results<-create_results_df(time)

#Read in single simulated data set
infile<-paste0("single_sim_observed_data.csv")
d<-read.csv(file=infile, sep = " ")

#Generate bootstrap samples
d.bs<-samp.bootstrap(nrow(d), B)

#Run bootstrap
for(j in 1:B){
  
  # Create bootstrap sample
  d_boot <- d[d.bs[,j], ]
  
  # Run estimators
  ice.bsresults[j,]<-as.numeric(ice_gcomp_run(d_boot, plan=plan_list))
  standard.bsresults[j,]<-as.numeric(standard_gcomp_run(d_boot, plan=plan_list, n_simulation=10000))     #n_simulation = number of iterations for monte carlo simulation   
      
  # CLEANUP: Remove bootstrap sample
  rm(d_boot)

  # Garbage collection every iterations
  gc(verbose = FALSE)
} 

# Loop through each method
for(estimator in estimators){

# Get the bsresults for this method
  bsresults <- get(paste0(estimator, ".bsresults"))
      
  new_row <- fill_results_row(bsresults)
      
  # Append to the appropriate results data frame
  results <- get(paste0(estimator, "_results"))
  assign(paste0(estimator, "_results"), rbind(results, new_row))
      
  # CLEANUP: Remove temporary variables
  rm(results, new_row)
}
  
# CLEANUP: Remove bootstrap results matrices
rm(ice.bsresults, standard.bsresults)
gc(verbose = FALSE)
  
ice_out<-paste0("single_sim_results_ice.csv")
standard_out<-paste0("single_sim_results_standard.csv")

write.table(ice_results,file=ice_out)
write.table(standard_results,standard_out)
