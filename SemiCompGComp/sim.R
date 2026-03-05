# Full simulation for recurrent and semi-competing event g-computation
# Alena Sorensen D'Alessio
# March 1, 2026


library(resample)
library(nnet)
library(stats)
library(betareg)
library(tidyverse)

setwd("YOUR DIRECTORY HERE")

source("multistate_DGM.R")
source("semicomp_gcomp_estimator_functions.R")
source("alternative_estimator_functions.R")
source("sim_functions.R")

#Designed to be run in parallel (suggest 100 parallel jobs)
args <- commandArgs(trailingOnly = TRUE) 
taskid <- as.numeric(args[1])
set.seed(taskid*5000)

#can ignore warnings as they are due to not having all contrasts of C(t) available during predict() because we haven't simulated through all time points yet 

#Simulation with bootstrap
n <- 2000 #number of observations (500, 2000)
B <- 500 #500 bootstraps
iterations <- 20 #20 iterations per task (with 100 parallel jobs, 2,000 iterations)

time<-3 #3 follow-up times after baseline
estimators<-c("alt1","alt2","ice","standard") #list of estimators
state_count<-3 #number of states (No hypertension, Hypertension, Death)
event<-c(1,2,3) #number corresponding to (no event, intermediate event, terminal event)
plan_list<-c("c(0,0,0)","c(1,1,1)") #untreated, treated plan
prefixes <- c("TrNon","TrInt", "TrTer", "NtrNon","NtrInt", "NtrTer")  #TrNon = treated, no intermediate or terminal event; TrInt = treated, intermediate event; TrTer = treated, terminal; NtrNon = Not treated, no intermediate or terminal event; NtrInt = Not treated, intermediate event; NtrTer = Not treated, terminal

#Create empty data frame for the bootstrap results
alt1.bsresults <- matrix(nrow = B, ncol = time*state_count*2) #Naive estimate ncol=time*state_count*2 (2 because we calculate mean and se for each measure)
alt2.bsresults <- matrix(nrow = B, ncol = time*state_count*2) #Naive estimate ncol=time*state_count*2 (2 because we calculate mean and se for each measure)
ice.bsresults <- matrix(nrow = B, ncol = time*state_count*2) # ICE
standard.bsresults <- matrix(nrow = B, ncol = time*state_count*2) # Monte-Carlo Flat

#Name columns of bootstrap results (not used but helps when checking mid-steps)
result_matrix_names <- unlist(lapply(prefixes, function(prefix) {
  unlist(lapply(1:time, function(i) {
    paste0(prefix, "_T", i) 
  }))
}))

dimnames(alt1.bsresults) <- list(NULL, result_matrix_names)
dimnames(alt2.bsresults) <- list(NULL, result_matrix_names)
dimnames(ice.bsresults) <- list(NULL, result_matrix_names)
dimnames(standard.bsresults) <- list(NULL, result_matrix_names)

#Create data frame for final results
alt1_results<-create_results_df(time)
alt2_results<-create_results_df(time)
ice_results<-create_results_df(time)
standard_results<-create_results_df(time)

#Run simulation with predefined number of iterations
for(i in 1:iterations){ 
  
  #Generate data with n observations
  d.potential<-dgm.potential(n) 
  d<-dgm.observed(d.potential) 
  
  #CLEANUP: Remove potential outcome data immediately after creating observed
  rm(d.potential)
  gc(verbose = FALSE)
  
  #Generate bootstrap samples
  d.bs<-samp.bootstrap(nrow(d), B)
  
  #Run bootstrap
  for(j in 1:B){
    # Create bootstrap sample
    d_boot <- d[d.bs[,j], ]
    
    #Run estimators
    alt1.bsresults[j,]<-as.numeric(alternative_estimator_run(d_boot,1, plan=plan_list))
    alt2.bsresults[j,]<-as.numeric(alternative_estimator_run(d_boot,2, plan=plan_list))
    ice.bsresults[j,]<-as.numeric(ice_gcomp_run(d_boot, plan=plan_list))
    standard.bsresults[j,]<-as.numeric(standard_gcomp_run(d_boot, plan=plan_list, n_simulation=10000)) #n_simulation = number of iterations for monte carlo simulation   
    
    #CLEANUP: Remove bootstrap sample
    rm(d_boot)
    
    #Garbage collection every iterations
    gc(verbose = FALSE)
  }
  
  #CLEANUP: Remove bootstrap indices and observed data
  rm(d.bs, d)

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
  
  # CLEANUP: Garbage collection after each simulation
  gc(verbose = FALSE)
}

# CLEANUP: Remove bootstrap results matrices
rm(alt1.bsresults, alt2.bsresults, ice.bsresults, standard.bsresults)
gc(verbose = FALSE)

#Save results
alt1_out    <-paste0("sim_output/sim_results_n", n, "_alt1_",taskid,".csv")
alt2_out    <-paste0("sim_output/sim_results_n", n, "_alt2_",taskid,".csv")
ice_out     <-paste0("sim_output/sim_results_n", n, "_ice_",taskid,".csv")
standard_out<-paste0("sim_output/sim_results_n", n, "_standard_",taskid,".csv")

write.table(alt1_results,file=alt1_out)
write.table(alt2_results,file=alt2_out)
write.table(ice_results,file=ice_out)
write.table(standard_results,standard_out)
