#####################################################################################################################
# Novel g-computation algorithms for time-varying actions with recurrent and semi-competing events
# D'Alessio et al.
# March 1, 2026

#Functions for proposed g-computation estimators
#####################################################################################################################


# Function to apply rmultinom to each row and handle missing
sample_multinomial <- function(prob_row) { #Currently not supported if all events are terminal
  if (is.na(prob_row[1])) return(NA)
  
  # Ensure prob_row has length 3 (for outcomes 1, 2, 3) (if less than 3 states, assume no terminal event occurred)
  if (length(prob_row) < 3) {
    full_prob <- rep(0, 3)
    names(full_prob) <- c("1", "2", "3")
    full_prob[1] <- 1-prob_row #set complement of intermediate event prob
    full_prob[2] <- prob_row #set intermediate event prob, terminal prob will remain zero
    prob_row <- full_prob
  }
  
  sample(1:3, size = 1, prob = prob_row)
}

# Function to apply rbinom to each row and handle missing
sample_binomial <- function(prob) {
  if(is.na(prob)){return(NA)}
  rbinom(1,1,prob)
}

get_colnames <- function(time){
  prefixes <- c("TrNon","TrInt", "TrTer", "NtrNon","NtrInt", "NtrTer")  #TrNon = treated, no intermediate or terminal event; TrInt = treated, intermediate event; TrTer = treated, terminal; NtrNon = Not treated, no intermediate or terminal event; NtrInt = Not treated, intermediate event; NtrTer = Not treated, terminal
  colnames <- unlist(lapply(prefixes, function(prefix) {
    unlist(lapply(1:time, function(i) {
      paste0(prefix, "_T", i) #format of TH_T1, TH_T2, etc.
    }))
  }))
  return(colnames)
}

#Standard g-computation estimator

standard_gcomp_simulation<-function(n_simulation,plan, dat,models){
  
  #get number of follow-ups
  plan<-eval(parse(text=plan))
  time<-length(plan)
  
  #set.seed(120225)  # Set seed for reproducibility
  
  #Sample with replacement *n_simulation* observations 
  dat.samp<- dat %>% sample_n(n_simulation, replace = TRUE)
  
    #Use L0 as start of sim data
    sim.dat<-data.frame(dat.samp)%>%select(L0)
    sim.dat$A0<-plan[[1]] 
    
    # CLEANUP: Remove dat.samp after extracting L0
    rm(dat.samp)
    
    #for all follow-ups not including the final outcome assessment
    for (t in 1:(time-1)) {
     A_var <- paste0("A", t)
     L_var <- paste0("L", t)
     Y_var <- paste0("Y", t)
     L_model <- paste0("L", t, ".model")
     Y_model <- paste0("Y", t, ".model")
     
     #predict Y at time t under plan 
     Yhat_pred<-predict(models[[Y_model]],sim.dat, type="prob")
     Yhat_pred<-as.matrix(Yhat_pred) #ensure Yhat_pred is a matrix (it is not in the case of no terminal events - sample_multinomial handles this case)
     Yhat <- apply(Yhat_pred, 1, sample_multinomial) #get outcome (1=no intermediate or terminal event, 1=intermediate event, 2=terminal event)
     sim.dat[[Y_var]]<-Yhat
     
     # CLEANUP: Remove prediction matrix
     rm(Yhat_pred, Yhat)
     
     #get indices for those still alive 
     alive <- sim.dat[[Y_var]] != event[3] & !is.na(sim.dat[[Y_var]])
     Lhat <- rep(NA, nrow(sim.dat))
     
     #predict L1 under plan 
     Lhat_pred<-predict.glm(models[[L_model]],sim.dat[alive,,drop=FALSE], type="response")
     Lhat[alive]<-rbinom(sum(alive),1,Lhat_pred)
     sim.dat<-sim.dat %>% mutate(!!L_var := case_when(.data[[Y_var]] == event[3] | is.na(.data[[Y_var]]) ~ NA,
                                                .default =Lhat),
                                 !!A_var := case_when(.data[[Y_var]] == event[3] | is.na(.data[[Y_var]]) ~ NA,
                                                .default =plan[[t+1]]))
     
     # CLEANUP
     rm(Lhat_pred, Lhat)
    }
    
    #predict Y at the final follow-up time under plan 
    Y_final_var <- paste0("Y", time)
    Y_final_model <- paste0("Y", time, ".model")
    
    Yhat_final_pred<-predict(models[[Y_final_model]],sim.dat, type="prob")
    Yhat_final_pred<-as.matrix(Yhat_final_pred) #ensure Yhat_pred is a matrix (it is not in the case of no terminal events - sample_multinomial handles this case)
    Yhat_final <- apply(Yhat_final_pred, 1, sample_multinomial) #get outcome (1=no intermediate or terminal event, 1=intermediate event, 2=terminal event)
    sim.dat[[Y_final_var]]<-Yhat_final
    
    # CLEANUP: Remove prediction matrix
    rm(Yhat_final_pred, Yhat_final)
    
    #Make composite outcome for each time point that reflects terminal event at any previous follow-up or uses the current Y
    for (t in 1:time) {
      Y_cols <- paste0("Y", 1:t) #get all Y's
      Y_current <-  paste0("Y", t)
      Y_fix <- paste0("Y", t, "_fix")
      if (t==1){sim.dat[[Y_fix]]=sim.dat[[Y_current]]}
      else{
        terminal <- paste0("terminal_by_Y", t)
        sim.dat[[terminal]] <- rowSums(sim.dat[, Y_cols] == event[3], na.rm=TRUE) > 0
        sim.dat[[Y_fix]]<-ifelse(sim.dat[[terminal]], event[3], sim.dat[[Y_current]]) #fix Y3 so that terminal event is indicated if individual had event at any time point
      }
    }
    
    #Initialize lists for intermediate and terminal proportion. Will be filled with values for each time point
    none_prop<-list()
    intermed_prop<-list()
    terminal_prop<-list()

    #Calculate proportion of intermediate and terminal events at each time point using the corrected outcome variables
    for (t in 1:time) {
      Y_var <- paste0("Y", t, "_fix")

      p.none    <-sum(sim.dat[[Y_var]]== event[1],na.rm=TRUE)/nrow(sim.dat) 
      p.intermed<-sum(sim.dat[[Y_var]]== event[2],na.rm=TRUE)/nrow(sim.dat) 
      p.terminal<-sum(sim.dat[[Y_var]]== event[3],na.rm=TRUE)/nrow(sim.dat) 
      
      none_prop[[t]]<-p.none
      intermed_prop[[t]]<-p.intermed
      terminal_prop[[t]]<- p.terminal
    }
    #Return proportion of Y @ each time as list of lists
    return(c(none_prop, intermed_prop, terminal_prop))
}

#run standard g-comp estimator
standard_gcomp_run<-function(dat, plan_list, n_simulation){

  # Model L1 given L0 and A0 
  L1.model <- glm(L1 ~ A0 + L0, data = dat, family = "binomial"(link = "logit"))
  # Model L1 given L0 and A0 
  L2.model <- glm(L2 ~ A1 + L1, data = dat, family = "binomial"(link = "logit"))
  # Model Y1 given L0 and A0 
  Y1.model <- multinom(Y1 ~ A0 + L0, data = dat, trace=FALSE)
  # Model Y2 given L1, L0, A1, and A0 
  Y2.model <- multinom(Y2 ~ A1 + L1 + Y1, data = dat, trace=FALSE)
  # Model Y3 given L1, L0, A1, and A0 
  Y3.model <- multinom(Y3 ~ A2 + L2 + Y2, data = dat, trace=FALSE)
  
  models <- list(Y1.model = Y1.model, Y2.model = Y2.model, Y3.model = Y3.model, 
                 L1.model=L1.model,L2.model=L2.model)
  
  #returns proportion of intermediate and terminal events at each time under treatment plan
  #order: intermediate.T1, intermediate.T2, intermediate.T3, terminal.T1, terminal.T2, terminal.T3
  treated<-standard_gcomp_simulation(n_simulation=n_simulation, plan=plan_list[2], dat=dat, models=models)   
  not_treated<-standard_gcomp_simulation(n_simulation=n_simulation, plan=plan_list[1], dat=dat, models=models) 
  
  # CLEANUP: Remove models
  rm(models, L1.model, L2.model, Y1.model, Y2.model, Y3.model)
  gc(verbose = FALSE)
  
  df<-data.frame(c(treated, not_treated))
  
  #add column names
  names(df)<-get_colnames(time)
  
  #Get treated and not treated proportions as one row data frame
  return(df) 
  
}
      
#ICE g-computation estimator

ice_gcomp_simulation<-function(dat, plan, model_formulas){
  
  #get number of follow-ups
  plan<-eval(parse(text=plan))
  time <- length(plan)
  
  #Restrict to those not censored at final time point
  C_final_var <- paste0("C", time)
  C_var <- paste0("C", time-1) #previous time
  A_var <- paste0("A", time-1)
  
  dat.sub<-dat%>%filter(.data[[C_final_var]]==0 | is.na(.data[[C_final_var]])) #Restrict to those uncensored at final time (can have terminal event)
  
  Y_final_model<-multinom(as.formula(paste0("Y", time, "~",model_formulas[time])),data=dat.sub) #Fit a regression model for Y at the final time using pre-specified formula for all observations not censored at end
  dat.sub <- dat %>% mutate(!!A_var:=plan[time])%>%filter(.data[[C_var]]==0| is.na(.data[[C_var]]))   #Set plan for previous time and restrict to those uncensored at previous time (can have terminal event)
  yhat <- as.matrix(predict(Y_final_model, newdata = dat.sub, type="prob")) #Generate predicted values of previous Y for the plan of interest (A at previous time is set to plan)
  
  #fix multinom output in case where there are no terminal events (if only one column results, that is the intermediate event prob, so add terminal event column with 0's and no event column with complement of intermediate event prob)
  if (ncol(yhat)<3) {
    missing_col<-matrix(0L, nrow = nrow(yhat), ncol = 1) 
    yhat <- cbind(missing_col, yhat, missing_col)
    yhat[,1] = 1-yhat[,2]
  }
  
  # CLEANUP
  rm(Y_final_model)
  
  #working backward from the (final time - 1) to 0 
  for (l in 1:(time-1)) { 
    t=time-l
    A_var <- paste0("A", t-1)
    C_var <- paste0("C", t-1)
    Y_var <- paste0("Y", t)
    
    #Replace predicted probability for those with terminal event at previous time with 100% probability of terminal event
    terminal_truth <- which(!is.na(dat.sub[[Y_var]]) & dat.sub[[Y_var]] == event[3])
    if(length(terminal_truth)>0){
      yhat[terminal_truth, ] <- matrix(c(0,0,1), 
                               nrow = length(terminal_truth), 
                               ncol = 3, 
                               byrow = TRUE) #gives warning if no one dies
    }
    assign(paste0("yhat", t), yhat)
    Y_model<-multinom(as.formula(paste0("yhat",t,"~",model_formulas[t])),data=dat.sub) #Fit a regression model for predicted Y at current time conditional on previous variables using predefined formula for all observations uncensored at current time
    
    if (t==1) dat.sub <- dat %>% mutate(!!A_var:=plan[t]) #no censoring at baseline
    if (t>1)  dat.sub <- dat %>% mutate(!!A_var:=plan[t])%>%filter(.data[[C_var]]==0| is.na(.data[[C_var]]))    #Set plan for previous time and restrict to those uncensored at previous time or current time if t=1
    yhat<-as.matrix(predict(Y_model,dat.sub, type="prob")) #Generate predicted values of Y1 for the plan of interest (e.g., A1=plan[2]) and uncensored at C1 
    
    #fix multinom output in case where there are no terminal events
    if (ncol(yhat)<3) {
      missing_col<-matrix(0L, nrow = nrow(yhat), ncol = 1) 
      yhat <- cbind(missing_col, yhat, missing_col)
      yhat[,1] = 1-yhat[,2]
    }
    
    # CLEANUP
    rm(Y_model)
  }

  assign(paste0("yhat", 0), yhat)
  rm(yhat)
  
  #No need to resolve probabilities at baseline time, just take means
  
  #Calculate proportion of intermediate and terminal events
  p.none<-mean(yhat0[,event[1]])
  p.intermed<-mean(yhat0[,event[2]])
  p.terminal<-mean(yhat0[,event[3]])
  
  result<-list(p.none, p.intermed, p.terminal)

  return(result)
}


#run ICE g-comp estimator
ice_gcomp_run<-function(dat, plan_list){
  
  rhs<-c("A0+L0","A1+L1+Y1", "A2+L2+Y2") # RHS of outcome equations. Order is t=1,2,3,...,etc.
  
  #Run simulation for each plan of interest
  treated<-ice_gcomp_simulation(dat, plan=plan_list[2], model_formulas=rhs)
  not_treated<-ice_gcomp_simulation(dat, plan=plan_list[1], model_formulas=rhs)
  
  #ICE does not have mid-point data so need to make sure to fill with NA's
  df<-data.frame(matrix(nrow=1, ncol=(time*state_count*2)))
  colnames(df)<-get_colnames(time)
  
  #fill with ICE results for final time point
  df[1,time*1]=treated[1]
  df[1,time*2]=treated[2]
  df[1,time*3]=treated[3]
  
  df[1,time*4]=not_treated[1]
  df[1,time*5]=not_treated[2]
  df[1,time*6]=not_treated[3]
  
  
  #Get treated and not treated proportions as one row dataframe
  return(df) 
}
