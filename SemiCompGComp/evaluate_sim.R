# Evaluation of full simulation for recurrent and semi-competing event g-computation
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


parallel=100 #100 parallel jobs
time=3 #3 time points with outcome after baseline
estimators<-c("alt1","alt2","ice","standard") #list of estimators
n=2000 #sample size
event<-c(1,2,3) #number corresponding to (intermediate event, terminal event)
state_list<-c("Non","Int","Ter")
plan_list<-c("c(0,0,0)","c(1,1,1)") #untreated, treated plan

#Import
alt1_results<-data.frame()
alt2_results<-data.frame()
ice_results<-data.frame()
standard_results<-data.frame()

for (estimator in estimators) {

  results_stack<-data.frame()
  
  for (i in 1:parallel) { 
    infile<-paste0("Output/Parallel/030326_multistate_g_comp_sim_results_n", n, "_",estimator,"_",i,".csv")
    results<-read.csv(file=infile, sep = " ")
    
    results_stack <- rbind(results_stack, results)
  }
  
  # Assign data frame back to the original name
  assign(paste0(estimator,"_results"), results_stack, envir = .GlobalEnv)
  rm(results_stack, results)
}

#get truth
n=1000000
set.seed(103125)
potential<-dgm.potential(n)
observed<-dgm.observed(potential) #error comes from missing censor data
truth<-get_truth(observed, plan_list)

# CLEANUP: Remove large truth truth dataset
rm(potential, observed)
gc(verbose = FALSE)

#calculate evaluation metrics

#CI Coverage
for (estimator in estimators) {
 #estimator<-"ice"
  dataset_name <- paste0(estimator, "_results")
  temp_df <- get(dataset_name)
  #t=3
  for (t in 1:time) {
    #add confidence intervals to each sim for CI coverage calculation
    # Get the data frame by its constructed name
    for(state in state_list){
    #state="Non"
    if((estimator == "alt1" | estimator == "alt2" | estimator == "ice" ) & t < time){ #alternative estimators and ice only have final time estimates so set other times to NA
      
      temp_df[[paste0("PD_", state, "_T", t, "_lci")]] <- NA
      temp_df[[paste0("PD_", state, "_T", t, "_uci")]] <- NA

      temp_df[[paste0("PD_", state, "_T", t, "_ci_coverage")]] <- NA

      temp_df[[paste0("PD_", state, "_T", t, "_bias")]] <- NA
    } else{
      # Add new columns
      
      if(estimator == "alt2" & state == "Ter"){ #alternative estimator 2 censors deaths so set terminal state values to NA
        
        temp_df[[paste0("PD_", state, "_T", t, "_lci")]] <- NA
        temp_df[[paste0("PD_", state, "_T", t, "_uci")]] <- NA
        
        temp_df[[paste0("PD_", state, "_T", t, "_ci_coverage")]] <- NA
        
        temp_df[[paste0("PD_", state, "_T", t, "_bias")]] <- NA
        
      }else{
      temp_df[[paste0("PD_", state, "_T", t, "_lci")]] <- temp_df[[paste0("PD_", state, "_T", t, "_mean")]]-1.96*temp_df[[paste0("PD_", state, "_T", t, "_se")]]
      temp_df[[paste0("PD_", state, "_T", t, "_uci")]] <- temp_df[[paste0("PD_", state, "_T", t, "_mean")]]+1.96*temp_df[[paste0("PD_", state, "_T", t, "_se")]]

      temp_df[[paste0("PD_", state, "_T", t, "_ci_coverage")]] <- ifelse(temp_df[[paste0("PD_", state, "_T", t, "_lci")]]<=(truth[[paste0("Tr", state, "_T", t)]]-truth[[paste0("Ntr", state, "_T", t)]]) & 
                                                                 temp_df[[paste0("PD_", state, "_T", t, "_uci")]]>=(truth[[paste0("Tr", state, "_T", t)]]-truth[[paste0("Ntr", state, "_T", t)]]),1,0)
    
      temp_df[[paste0("PD_", state, "_T", t, "_bias")]] <- temp_df[[paste0("PD_", state, "_T", t, "_mean")]]-(truth[[paste0("Tr", state, "_T", t)]]-truth[[paste0("Ntr", state, "_T", t)]])
      }
    }
    # Assign data frame back to the original name
    assign(dataset_name, temp_df, envir = .GlobalEnv)
  }
  }
}

estimates <- list()
bias <- list()
ese <- list()
avg.se <- list()
mse <- list()
avg.lci <- list()
avg.uci <- list()
cld <- list()
coverage <- list()

for (estimator in estimators) {

  for (t in 1:time) {
   
    # Get the data frame by its constructed name
    for(state in state_list){
      if(estimator=="ice" & t < time){ #ice only has final time estimates so set other times to NA

      estimates[[paste0("PD.", state, ".T", t, ".", estimator)]] <- NA 

      bias[[paste0("bias.PD.", state, ".T", t, ".", estimator)]] <-NA

      ese[[paste0("ese.PD.", state, ".T", t, ".", estimator)]] <-NA

      avg.se[[paste0("avg.se.PD.", state, ".T", t, ".", estimator)]] <-NA

      mse[[paste0("mse.PD.", state, ".T", t, ".", estimator)]] <-NA

      avg.lci[[paste0("avg.lci.PD.", state, ".T", t, ".", estimator)]] <-NA

      avg.uci[[paste0("avg.uci.PD.", state, ".T", t, ".", estimator)]] <-NA

      cld[[paste0("cld.PD.", state, ".T", t, ".", estimator)]] <-NA

      coverage[[paste0("coverage.PD.", state, ".T", t, ".", estimator)]] <-NA

      } else{
      
      results <- get(paste0(estimator, "_results"))
      
      # Add new columns
      estimates[[paste0("Tr", state, ".T", t, ".", estimator)]] <- 
        mean(results[[paste0("Tr", state, "_T", t, "_mean")]]) #Prevalence of hypertension among treated at time t

      estimates[[paste0("Ntr", state, ".T", t, ".", estimator)]] <- 
        mean(results[[paste0("Ntr", state, "_T", t, "_mean")]]) #Prevalence of hypertension among untreated at time t

      estimates[[paste0("PD.", state, ".T", t, ".", estimator)]] <- 
        mean(results[[paste0("PD_", state, "_T", t, "_mean")]]) #Prevalence difference of hypertension at time t

      bias[[paste0("bias.PD.", state, ".T", t, ".", estimator)]] <- mean(results[[paste0("PD_", state, "_T", t, "_bias")]])

      ese[[paste0("ese.PD.", state, ".T", t, ".", estimator)]] <- 
        sd(results[[paste0("PD_", state, "_T", t, "_mean")]])

      avg.se[[paste0("avg.se.PD.", state, ".T", t, ".", estimator)]] <- 
        mean(results[[paste0("PD_", state, "_T", t, "_se")]])

      mse[[paste0("mse.PD.", state, ".T", t, ".", estimator)]] <- 
        bias[[paste0("bias.PD.", state, ".T", t, ".", estimator)]]^2+ese[[paste0("ese.PD.", state, ".T", t, ".", estimator)]]^2

      avg.lci[[paste0("avg.lci.PD.", state, ".T", t, ".", estimator)]] <- 
        estimates[[paste0("PD.", state, ".T", t, ".", estimator)]]-1.96*avg.se[[paste0("avg.se.PD.", state, ".T", t, ".", estimator)]]
            
      avg.uci[[paste0("avg.uci.PD.", state, ".T", t, ".", estimator)]] <- 
        estimates[[paste0("PD.", state, ".T", t, ".", estimator)]]+1.96*avg.se[[paste0("avg.se.PD.", state, ".T", t, ".", estimator)]]

      cld[[paste0("cld.PD.", state, ".T", t, ".", estimator)]] <- 
        avg.uci[[paste0("avg.uci.PD.", state, ".T", t, ".", estimator)]] - avg.lci[[paste0("avg.lci.PD.", state, ".T", t, ".", estimator)]]
            
      coverage[[paste0("coverage.PD.", state, ".T", t, ".", estimator)]] <- 
        mean(results[[paste0("PD_", state, "_T", t, "_ci_coverage")]])
 
            
      # CLEANUP: Remove results variable
      rm(results)
      }
      }
  }
}

#get PD estimates at each time point
Estimator_name<-c("Truth","Alternative Estimator 1","Alternative Estimator 2", "ICE G-comp", "Standard G-comp")
for (t in 1:time) {
  # Store data frames for each state at this time point
  state_dfs <- list()
  
  for(state in state_list){
  PD_Estimate<-c(truth[[paste0("Tr", state, "_T", t)]]-truth[[paste0("Ntr", state, "_T", t)]],estimates[[paste0("PD.", state, ".T", t, ".alt1")]],estimates[[paste0("PD.", state, ".T", t, ".alt2")]],estimates[[paste0("PD.", state, ".T", t, ".ice")]],estimates[[paste0("PD.", state, ".T", t, ".standard")]])
  PD_Lower_CI <- c(NA, avg.lci[[paste0("avg.lci.PD.", state, ".T", t, ".alt1")]], avg.lci[[paste0("avg.lci.PD.", state, ".T", t, ".alt2")]],avg.lci[[paste0("avg.lci.PD.", state, ".T", t, ".ice")]],avg.lci[[paste0("avg.lci.PD.", state, ".T", t, ".standard")]])
  PD_Upper_CI <- c(NA, avg.uci[[paste0("avg.uci.PD.", state, ".T", t, ".alt1")]],avg.uci[[paste0("avg.uci.PD.", state, ".T", t, ".alt2")]],avg.uci[[paste0("avg.uci.PD.", state, ".T", t, ".ice")]],avg.uci[[paste0("avg.uci.PD.", state, ".T", t, ".standard")]])
  PD_CLD <- c(NA, cld[[paste0("cld.PD.", state, ".T", t, ".alt1")]],cld[[paste0("cld.PD.", state, ".T", t, ".alt2")]],cld[[paste0("cld.PD.", state, ".T", t, ".ice")]],cld[[paste0("cld.PD.", state, ".T", t, ".standard")]])
  PD_Coverage <- c(NA, coverage[[paste0("coverage.PD.", state, ".T", t, ".alt1")]],coverage[[paste0("coverage.PD.", state, ".T", t, ".alt2")]],coverage[[paste0("coverage.PD.", state, ".T", t, ".ice")]],coverage[[paste0("coverage.PD.", state, ".T", t, ".standard")]])
  PD_Bias<-c(NA, bias[[paste0("bias.PD.", state, ".T", t, ".alt1")]],bias[[paste0("bias.PD.", state, ".T", t, ".alt2")]],bias[[paste0("bias.PD.", state, ".T", t, ".ice")]],bias[[paste0("bias.PD.", state, ".T", t, ".standard")]])
  PD_ESE<-c(NA, ese[[paste0("ese.PD.", state, ".T", t, ".alt1")]],ese[[paste0("ese.PD.", state, ".T", t, ".alt2")]],ese[[paste0("ese.PD.", state, ".T", t, ".ice")]],ese[[paste0("ese.PD.", state, ".T", t, ".standard")]])
  PD_MSE<-c(NA, mse[[paste0("mse.PD.", state, ".T", t, ".alt1")]],mse[[paste0("mse.PD.", state, ".T", t, ".alt2")]],mse[[paste0("mse.PD.", state, ".T", t, ".ice")]],mse[[paste0("mse.PD.", state, ".T", t, ".standard")]])
  PD_ASE<-c(NA, avg.se[[paste0("avg.se.PD.", state, ".T", t, ".alt1")]],avg.se[[paste0("avg.se.PD.", state, ".T", t, ".alt2")]],avg.se[[paste0("avg.se.PD.", state, ".T", t, ".ice")]],avg.se[[paste0("avg.se.PD.", state, ".T", t, ".standard")]])
  PD_SER<-c(NA, avg.se[[paste0("avg.se.PD.", state, ".T", t, ".alt1")]]/ese[[paste0("ese.PD.", state, ".T", t, ".alt1")]],avg.se[[paste0("avg.se.PD.", state, ".T", t, ".alt2")]]/ese[[paste0("ese.PD.", state, ".T", t, ".alt2")]], avg.se[[paste0("avg.se.PD.", state, ".T", t, ".ice")]]/ese[[paste0("ese.PD.", state, ".T", t, ".ice")]], avg.se[[paste0("avg.se.PD.", state, ".T", t, ".standard")]]/ese[[paste0("ese.PD.", state, ".T", t, ".standard")]])
  
  temp_df<-data.frame(Estimator_name, PD_Estimate, PD_Lower_CI, PD_Upper_CI, PD_CLD, PD_Coverage, PD_Bias, PD_ESE, PD_MSE, PD_ASE, PD_SER)
  
  temp_df <- temp_df %>%
    pivot_longer(-1, names_to = "variable", values_to = "value") %>%
    pivot_wider(names_from = 1, values_from = value)
  
  temp_df$variable = paste0(state,"_", temp_df$variable)
  
  # Store in list instead of immediately assigning
  state_dfs[[state]] <- temp_df
  }
  
  # Combine all states for this time point
  combined_df <- do.call(rbind, state_dfs)
  
  # Assign the combined data frame
  assign(paste0("estimator_comparison_T", t), combined_df, envir = .GlobalEnv)

}

t1_out<-paste0("Output/Parallel/030336_multistate_g_comp_sim_results_n", n, "_T1.csv")
t2_out<-paste0("Output/Parallel/030336_multistate_g_comp_sim_results_n", n, "_T2.csv")
t3_out<-paste0("Output/Parallel/030336_multistate_g_comp_sim_results_n", n, "_T3.csv")

write.table(estimator_comparison_T1,file=t1_out)
write.table(estimator_comparison_T2,file=t2_out)
write.table(estimator_comparison_T3,file=t3_out)
