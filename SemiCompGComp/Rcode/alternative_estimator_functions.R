#####################################################################################################################
# Novel g-computation algorithms for time-varying actions with recurrent and semi-competing events
# D'Alessio et al.
# March 1, 2026

#Functions for alternate estimators
#####################################################################################################################


#Alternate estimator 1
#Multistate g-computation estimator that only considers the action at baseline (i.e., an intent to treat estimator for a per protocol parameter). 
alternative_estimator_1<-function(dat, plan){

  #get number of follow-ups
  plan<-eval(parse(text=plan))
  time <- length(plan)
  
  C_final_var <- paste0("C", time)
  Y_final_var <- paste0("Y", time)
  A_var <- paste0("A", 0) #ITT so only assess baseline A
  rhs<-"A0+L0" #No accounting for time-varying so only uses baseline data
  
  #Update outcome at final time to reflect terminal events at prior time
  #Can identify terminal event occurrence if Y and C are missing at final time
  dat<-dat %>% mutate(!!Y_final_var := case_when(is.na(.data[[Y_final_var]]) & is.na(.data[[C_final_var]]) ~ event[3],
                                          .default = .data[[Y_final_var]]))
  
  #Restrict to those with an event status at final time point (have terminal event or are not censored) and on plan at baseline
  dat.sub<-dat%>%filter(!is.na(.data[[Y_final_var]]) & .data[[A_var]]==plan[1]) 
  
  #Fit a regression model for Y at the final time for all observations with outcome at final time
  Y_final_model<-multinom(as.formula(paste0("Y", time, "~",rhs)),data=dat.sub, family="binomial") 
  
  Y_final_model<-multinom(as.formula(paste0("Y", time, "~",rhs)),data=dat.sub) #Fit a regression model for Y at the final time using pre-specified formula for all observations not censored at end
  yhat <- as.matrix(predict(Y_final_model, newdata = dat.sub, type="prob")) #Generate predicted values of previous Y for the plan of interest (A at previous time is set to plan)
  
  #Initialize lists for intermediate and terminal proportion. Will be filled with values for each time point
  none_prop<-list()
  intermed_prop<-list()
  terminal_prop<-list()
  
  #Calculate proportion of intermediate and terminal events final point
  for (t in 1:time) {
    none_prop[[t]]<- NA 
    intermed_prop[[t]]<- NA 
    terminal_prop[[t]]<- NA 
    
    if (t==time){
      none_prop[[t]]<-mean(yhat[,event[1]])
      intermed_prop[[t]]<-mean(yhat[,event[2]])
      terminal_prop[[t]]<-mean(yhat[,event[3]])
    }
  }
  
  #Return proportion of Y at final time
  return(c(none_prop, intermed_prop, terminal_prop))
}

#Alternate estimator 2
#ICE g-computation estimator that accounts for time-varying confounding but censors the terminal event.
alternative_estimator_2<-function(dat, plan){
  
  rhs<-c("A0+L0", "A1+L1+Y1", "A2+L2+Y2") # RHS of outcome equations. Order is t=1,2,3,...,etc.
  
  #get number of follow-ups
  plan<-eval(parse(text=plan))
  time <- length(plan)
  
  #Update those who had terminal event to censored at each time. 
  for (t in 1:time) { 
    C_var <- paste0("C", t) 
    Y_var <- paste0("Y", t)
    
    dat<-dat %>% mutate(!!C_var := case_when(.data[[Y_var]] == event[3] | is.na(.data[[C_var]]) ~ 1,
                                             .default = .data[[C_var]]),
                        !!Y_var :=case_when(.data[[Y_var]] == event[3] ~ NA,
                                            .default = .data[[Y_var]]))
  }
  
  #subset sample to those uncensored at T1 (if censored at T1 that means no outcomes were ever observed and shouldn't be included in analysis)
  dat <- dat %>% filter(C1==0) 
  
  #Restrict to those not censored at final time point (terminal event will be treated as censoring)
  C_final_var <- paste0("C", time)
  Y_final_var <- paste0("Y", time)
  C_var <- paste0("C", time-1) #next time
  A_var <- paste0("A", time-1)
  Y_var <- paste0("Y", time-1)
  
  #Restrict to those uncensored at final time 
  dat.sub<-dat%>%filter(.data[[C_final_var]]==0) %>% mutate(!!Y_final_var := as.factor(.data[[Y_final_var]]))   #Update Y to a factor excluding the terminal event
  
  #Fit a regression model for Y at the final time using pre-specified formula for all observations not censored at end (terminal event is treated as censoring so logistic model needed)
  Y_final_model<-glm(as.formula(paste0("Y", time, "~",rhs[time])),data=dat.sub, family="binomial") 
  
  #Set plan for previous time and restrict to those uncensored at previous time
  dat.sub <- dat %>% mutate(!!A_var:=plan[time])%>%filter(.data[[C_var]]==0)
  
  #Generate predicted values of previous Y for the plan of interest (A at previous time is set to plan)
  yhat <- predict(Y_final_model, newdata = dat.sub, type="response") 
  
  # CLEANUP
  rm(Y_final_model)
  
  #working backward from the (final time - 1) to 0 
  for (l in 1:(time-1)) { 
    #l=2
    t=time-l
    A_var <- paste0("A", t-1)
    C_var <- paste0("C", t-1)
    Y_var <- paste0("Y", t)
    
    #Update Y to a factor excluding the terminal event.  Previously subset data to those uncensored and without terminal event at current time
    dat.sub <- dat.sub %>% mutate(!!Y_var := as.factor(.data[[Y_var]]))   
    
    #Fit a regression model for predicted Y at current time conditional on previous variables using predefined formula for all observations uncensored at current time
    assign(paste0("yhat", t), yhat)
    Y_model<-betareg(as.formula(paste0("yhat",t,"~",rhs[t])),data=dat.sub) 
    
    #Set plan for previous time
    if (t==1) dat.sub <- dat %>% mutate(!!A_var:=plan[t]) #No censoring at baseline
    if (t>1)  dat.sub <- dat %>% mutate(!!A_var:=plan[t])%>%filter(.data[[C_var]]==0)    #Restrict to those uncensored at previous time
    
    #Generate predicted values of previous Y for the plan of interest 
    yhat<-predict(Y_model,dat.sub) 
    
    # CLEANUP
    rm(Y_model)
  }
  
  assign(paste0("yhat", 0), yhat)
  rm(yhat)
  
  #Initialize lists for intermediate and terminal proportion. Will be filled with values for each time point
  none_prop<-list()
  intermed_prop<-list()
  terminal_prop<-list()
  
  #Calculate proportion of intermediate and terminal events at each time point using the corrected outcome variables
  for (t in 1:time) {
    none_prop[[t]]<- NA #ICE only provides estimate at final time
    intermed_prop[[t]]<- NA #ICE only provides estimate at final time
    terminal_prop[[t]]<- NA #terminal event treated as censoring
    
    if (t==time){
      none_prop[[t]]<-1-mean(yhat0) 
      intermed_prop[[t]]<-mean(yhat0) 
    }
  }
  #Return proportion of Y @ each time as list of lists
  return(c(none_prop, intermed_prop, terminal_prop))
}

#Run alternative estimators 
#version options = 1 or 2
alternative_estimator_run<-function(dat, version, plan_list){
  #Run estimator for each plan of interest
  if(version == 1){
    treated<-alternative_estimator_1(dat, plan=plan_list[2])
    not_treated<-alternative_estimator_1(dat, plan=plan_list[1])
  }else if(version == 2){
    treated<-alternative_estimator_2(dat, plan=plan_list[2])
    not_treated<-alternative_estimator_2(dat, plan=plan_list[1])
  }else{
    stop("Invalid version specified. Must be one of: '1' or '2'")
  }
  
  df<-data.frame(c(treated, not_treated))
  
  #add column names
  names(df)<-get_colnames(time)
  #Get treated and not treated proportions as one row dataframe
  return(df) 
}
