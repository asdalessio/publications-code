# Data generating mechanism for recurrent and semi-competing events
# Alena Sorensen D'Alessio
# March 1, 2026


##Helper functions

#Expit function
expit<-function(x){
  exp(x)/(1+exp(x))
}

#Calculate outcome probabilities
calculate_probs <- function(eq_vals) {
  denominator <- 1 + exp(eq_vals[1]) + exp(eq_vals[2])
  return(c(1/denominator, exp(eq_vals[1])/denominator, exp(eq_vals[2])/denominator))
}


#Outcome equations

#time=1
outcome_eq_t1<-function(a_0, L_0){
  hypertensive<-0.05-0.6*a_0-2*L_0
  dead<--1.75-0.6*a_0-2*L_0
  return(c(hypertensive, dead)) #returns list of hypertensive prob, dead prob
}

#time=2
outcome_eq_t2<-function(a_1, L_1, Y_1){
  hypertensive<- 0.1-0.8*a_1-2.2*L_1+0.5*Y_1
  dead       <--2-0.6*a_1-2*L_1+0.4*Y_1
  return(c(hypertensive, dead)) #returns list of hypertensive prob, dead prob
}

#time=3
outcome_eq_t3<-function(a_2, L_2, Y_2){
  hypertensive<- 0.3-0.9*a_2-2.2*L_2+0.5*Y_2 
  dead       <--2-0.6*a_2-2*L_2+0.4*Y_2
  return(c(hypertensive, dead)) #returns list of hypertensive prob, dead prob
}

#Get outcome function

#time=1
multinom_outcome_t1 <- function(a_0, L_0) {
  #Calculate values from outcome eq fucntions
  outcome_eq_result <- outcome_eq_t1(a_0, L_0) #first column is value for hypertensive, second column is dead
  
  # Calculate probabilities 
  probs <- calculate_probs(outcome_eq_result)
  
  #Get single, resolved outcome
  return(which(rmultinom(1,1,probs)==1))
}

#time=2
multinom_outcome_t2 <- function(a_1, L_1, Y_1){
  #Calculate values from outcome eq fucntions
  outcome_eq_result <- outcome_eq_t2(a_1, L_1, Y_1) #first column is value for hypertensive, second column is dead
  
  # Calculate probabilities 
  probs <- calculate_probs(outcome_eq_result)
  
  #Get single, resolved outcome
  return(which(rmultinom(1,1,probs)==1))
}

#time=3
multinom_outcome_t3 <- function(a_2, L_2, Y_2){
  #Calculate values from dgm outcome eq
  outcome_eq_result <- outcome_eq_t3(a_2, L_2, Y_2) #first column is value for hypertensive, second column is dead
  
  # Calculate probabilities 
  probs <- calculate_probs(outcome_eq_result)
  
  #Get single, resolved outcome
  return(which(rmultinom(1,1,probs)==1))
}

#Calculate potential outcomes

dgm.potential <-function(n){
  dat<-data.frame(
                  #L potential outcomes
                  #Bernoulli(0.5)
                  L0=rbinom(n,1,0.5))
  dat<-dat%>%mutate(
                  #L potential outcomes
                  #Bernoulli {expit(−1 − a0 + Li,0)}
                  L1.a0=rbinom(n,1,expit(-1-0+L0)), #(a0==0) 
                  L1.a1=rbinom(n,1,expit(-1-1+L0)), #(a0==1) 
                  
                  #Bernoulli {expit(−1 − a1 + Li,1(a0))}
                  L2.a0a0=rbinom(n,1,expit(-1-0+L1.a0)),  #(a0==0, a1==0)
                  L2.a0a1=rbinom(n,1,expit(-1-1+L1.a0)),  #(a0==0, a1==1) 
                  L2.a1a0=rbinom(n,1,expit(-1-0+L1.a1)),  #(a0==1, a1==0) 
                  L2.a1a1=rbinom(n,1,expit(-1-1+L1.a1)))  #(a0==1, a1==1)
  
  # Free up memory after L generation
  gc() 
                  
    #Y potential outcomes 
    #Y @ t=1
    dat$Y1.a0 = vapply(1:n, function(i) multinom_outcome_t1(0,dat$L0[i]), numeric(1))
    dat$Y1.a1 = vapply(1:n, function(i) multinom_outcome_t1(1,dat$L0[i]), numeric(1))
    
    gc()
    
    #Y @ t=2 
    dat$Y2.a0a0 = ifelse(dat$Y1.a0==3, 3, vapply(1:n, function(i) multinom_outcome_t2(0, dat$L1.a0[i], dat$Y1.a0[i]), numeric(1)))
    dat$Y2.a1a0 = ifelse(dat$Y1.a1==3, 3, vapply(1:n, function(i) multinom_outcome_t2(0, dat$L1.a1[i], dat$Y1.a1[i]), numeric(1)))
    dat$Y2.a0a1 = ifelse(dat$Y1.a0==3, 3, vapply(1:n, function(i) multinom_outcome_t2(1, dat$L1.a0[i], dat$Y1.a0[i]), numeric(1)))
    dat$Y2.a1a1 = ifelse(dat$Y1.a1==3, 3, vapply(1:n, function(i) multinom_outcome_t2(1, dat$L1.a1[i], dat$Y1.a1[i]), numeric(1)))
    
    gc()
    
    #Y @ t=3
    dat$Y3.a0a0a0 = ifelse(dat$Y2.a0a0==3, 3, vapply(1:n, function(i) multinom_outcome_t3(0, dat$L2.a0a0[i], dat$Y2.a0a0[i]), numeric(1)))
    dat$Y3.a1a0a0 = ifelse(dat$Y2.a1a0==3, 3, vapply(1:n, function(i) multinom_outcome_t3(0, dat$L2.a1a0[i], dat$Y2.a1a0[i]), numeric(1)))
    dat$Y3.a0a1a0 = ifelse(dat$Y2.a0a1==3, 3, vapply(1:n, function(i) multinom_outcome_t3(0, dat$L2.a0a1[i], dat$Y2.a0a1[i]), numeric(1)))
    dat$Y3.a0a0a1 = ifelse(dat$Y2.a0a0==3, 3, vapply(1:n, function(i) multinom_outcome_t3(1, dat$L2.a0a0[i], dat$Y2.a0a0[i]), numeric(1)))
    dat$Y3.a1a1a0 = ifelse(dat$Y2.a1a1==3, 3, vapply(1:n, function(i) multinom_outcome_t3(0, dat$L2.a1a1[i], dat$Y2.a1a1[i]), numeric(1)))
    dat$Y3.a0a1a1 = ifelse(dat$Y2.a0a1==3, 3, vapply(1:n, function(i) multinom_outcome_t3(1, dat$L2.a0a1[i], dat$Y2.a0a1[i]), numeric(1)))
    dat$Y3.a1a0a1 = ifelse(dat$Y2.a1a0==3, 3, vapply(1:n, function(i) multinom_outcome_t3(1, dat$L2.a1a0[i], dat$Y2.a1a0[i]), numeric(1)))
    dat$Y3.a1a1a1 = ifelse(dat$Y2.a1a1==3, 3, vapply(1:n, function(i) multinom_outcome_t3(1, dat$L2.a1a1[i], dat$Y2.a1a1[i]), numeric(1)))
    
    gc()
    
  return(dat)
}

#Get observed data

dgm.observed <-function(dat){
  n<-nrow(dat)
  dat<-dat %>% mutate(
                  #Observed treatment and confounder values
                  A0=rbinom(n,1,expit(1-2*L0)),  #Bernoulli {expit(1 − 2Li,0)}    
                  L1=case_when(A0==1 ~ L1.a1, .default=L1.a0),
                  A1=rbinom( n,1,expit(-1-L1+1.75*A0)), #Bernoulli {expit(−1 − Li,1 + 1.75Ai,0)}
                  L2=case_when(A0==1 & A1==1 ~ L2.a1a1,
                               A0==1 & A1==0 ~ L2.a1a0,
                               A0==0 & A1==1 ~ L2.a0a1,
                               .default=L2.a0a0),
                  A2=rbinom(n,1,expit(-1-L2+1.75*A1)), #Bernoulli {expit(−1 − Li,2 + 1.75Ai,1)}
                  
                  #Observed outcome values
                  Y1=case_when(A0==1 ~ Y1.a1,
                               .default=Y1.a0),
                  Y2=case_when(A0==1 & A1==1 ~ Y2.a1a1,
                               A0==1 & A1==0 ~ Y2.a1a0,
                               A0==0 & A1==1 ~ Y2.a0a1,
                               .default=Y2.a0a0),
                  Y3=case_when(A0==1 & A1==1 & A2==1 ~ Y3.a1a1a1,
                               A0==0 & A1==1 & A2==1 ~ Y3.a0a1a1,
                               A0==1 & A1==0 & A2==1 ~ Y3.a1a0a1,
                               A0==1 & A1==1 & A2==0 ~ Y3.a1a1a0,
                               A0==0 & A1==0 & A2==1 ~ Y3.a0a0a1,
                               A0==0 & A1==1 & A2==0 ~ Y3.a0a1a0,
                               A0==1 & A1==0 & A2==0 ~ Y3.a1a0a0,
                               .default=Y3.a0a0a0))
                  
  #Censoring
  dat<-dat %>% mutate(
                 C1=rbinom(n=n,size=1,prob=expit(-3-0.5*A0)), #C(i,1) ∼ Bernoulli {expit(−3 − 0.5Ai,0)} 
                 Y1=case_when(C1==1 ~ NA,
                              .default=Y1),
                 A1=case_when(C1==1 ~ NA,
                              !is.na(Y1) & Y1==3 ~ NA,
                              .default=A1),
                 L1=case_when(C1==1 ~ NA,
                               !is.na(Y1) & Y1==3 ~ NA,
                               .default=L1))%>%
                mutate(
             # Generate C2 only for those non-censored or dead at t=1
                  C2=case_when(C1==1 ~ 1,
                               Y1==3 ~ NA, 
                               .default = rbinom(n = n, size = 1, prob = expit(-3 - 0.5*A1))),
                  Y2=case_when(C2==1 ~ NA,
                               Y1==3 ~ NA,
                               .default = Y2),
                  A2=case_when(C1==1 ~ NA,
                               is.na(Y2) ~ NA,
                               !is.na(Y2) & Y2==3 ~ NA,
                               .default=A2),
                  L2=case_when(C1==1 ~ NA,
                               is.na(Y2) ~ NA,
                               !is.na(Y2) & Y2==3 ~ NA,
                               .default=L2)) %>%
                mutate(
              # Generate C3 only for those non-censored or dead at t=1 or t=2
                  C3=case_when(C1==1 | C2==1 ~ 1,
                               Y1==3 | Y2==3 ~ NA, 
                                .default = rbinom(n = n, size = 1, prob = expit(-3-0.5*A2))),
                  Y3=case_when(C3==1 ~ NA,
                               Y1==3 | Y2==3 ~ NA, 
                               .default = Y3)) 
  return(dat)
}
