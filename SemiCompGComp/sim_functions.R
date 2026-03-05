# Functions for simulation for recurrent and semi-competing event g-computation
# Alena Sorensen D'Alessio
# March 1, 2026


#Helper function to identify true proportions
truth_helper<-function(dat, plan){
  
  #get number of follow-ups
  plan<-eval(parse(text=plan))
  time <- length(plan)
  
  #Initialize lists for intermediate and terminal events proportion under each treatment plan. Will be filled with values for each time point
  none_prop<-list()
  intermed_prop<-list()
  terminal_prop<-list()
  
  for (t in 1:time) {
    # Generate the potential outcome pattern: "a1" or "a0" repeated 'time' number of times
    po.suffix <-paste0("a", plan[1:t], collapse = "") 
    
    p.none<-sum(dat[[paste0("Y", t, ".", po.suffix)]]==event[1],na.rm=TRUE)/nrow(dat) 
    p.intermed<-sum(dat[[paste0("Y", t, ".", po.suffix)]]==event[2],na.rm=TRUE)/nrow(dat) 
    p.terminal<-sum(dat[[paste0("Y", t, ".", po.suffix)]]==event[3],na.rm=TRUE)/nrow(dat) 
    
    none_prop[[t]]<-p.none
    intermed_prop[[t]]<-p.intermed
    terminal_prop[[t]]<- p.terminal
  }
  
  #Return proportions of Y @ each time as list of lists
  return(c(none_prop, intermed_prop, terminal_prop))
}

#Function to identify true proportions
get_truth<-function(dat, plan_list){
  
  #Get true proportions for each plan of interest
  treated<-truth_helper(dat, plan=plan_list[2])
  not_treated<-truth_helper(dat, plan=plan_list[1])
  
  df<-data.frame(c(treated, not_treated))
  
  #add column names
  names(df)<-get_colnames(time)
  
  #Get treated and not treated proportions as one row data frame
  return(df) 
}

#Simulation management functions

#Create the results data frame
create_results_df <- function(time) {
  all_measures<-get_result_colnames(time)
  # Create empty data frame with double type columns
  df <- setNames(
    data.frame(matrix(numeric(0),ncol = length(all_measures))),
    all_measures
  )
  return(df)
}

#Condense bootstrapped results into single row at the appropriate location
fill_results_row <- function(mat){
  col=1
  estimate_row=time*state_count*2
  temp<-matrix(nrow=1, ncol=(estimate_row*3))
  dimnames(temp) <- list(NULL, get_result_colnames(time))
  
  for (k in 1:estimate_row){ #number of rows per treated, not treated, prev difference groups
    temp[,col] <- mean(mat[,k]) 
    temp[,col+1] <- sd(mat[,k])
    col=col+2
  }
  
  contrast_row=estimate_row/2 #divide by number of contrasts
  for (k in 1:contrast_row){ #number of rows per treated, not treated, prev difference groups
    temp[,col] <- mean(mat[,k]-mat[,(k+state_count*time)])
    temp[,col+1] <- sd(mat[,k]-mat[,(k+state_count*time)])
    col=col+2
  }
  
  return(temp)
}

#Apply the appropriate column name for the simulation results data frame
get_result_colnames <-function(time){
  # Define prefixes and suffixes
  prefixes <- c("TrNon","TrInt", "TrTer", "NtrNon","NtrInt", "NtrTer")  #TrNon = treated, no intermediate or terminal event; TrInt = treated, intermediate event; TrTer = treated, terminal; NtrNon = Not treated, no intermediate or terminal event; NtrInt = Not treated, intermediate event; NtrTer = Not treated, terminal
  measures <- c("mean", "se")
  
  # Create column names for the first follow-up time 
  proportion <- unlist(lapply(prefixes, function(prefix) {
    unlist(lapply(1:time, function(i) {
      paste0(prefix, "_T", i, "_", measures) #format of TH_T1_mean, NTD_T1_mean, etc.
    }))
  }))
  
  # Create column names for proportion difference columns
  pd_prefixes <- c("PD_Non","PD_Int", "PD_Ter") #_Non = no intermediate or terminal event; _Int = intermediate event; _Ter = terminal
  prev_difference <- unlist(lapply(pd_prefixes, function(pd_prefix) {
    unlist(lapply(1:time, function(i) {
      paste0(pd_prefix, "_T", i, "_", measures)
    }))
  }))
  
  # Combine all column names
  all_measures <- c(proportion, prev_difference)
  return(all_measures)
}
