######################################################################
### This code provides all the required baseline chracteristics 
### and Specification of the distribution of baseline characteristics
### It includes 4 parameters:
### p_HIVMSM, p_IDU_start, p_IDU_stop, P_HIV
### Also it include the proportion of patients in each fibrosis stage (start)
### We calculate time to get HIV, which depends on the probability of getting HIV after HCV, and it can be a number between 3 month to maxTime, or it never happens. 
### We calculate the time of starting and stopping IDU, which depend on the probability of starting and stoping IDU respectively.
### We also have two function:
# - sr.fun to calculate the spontaneous recovery time
# - pw.eval.ext to calculate the background mortality
#####################################################################
#setwd("../Required_Data")
bl_number = 13 # number of baseline characteristics
statesNum = 52 # number of states
#proportion of patients in each fibrosis stages (F0, F1, F2, F3, F4) at simulation start (note: simulation starts at time of HCV infection)
start <- sample(c(1,2,3,4,5), cohortSize, c(0.8590, 0.1509, 0.0184, 0.00239, 0.00015), replace=T) 
p_HIVMSM = 0.5
p_IDU_start = 0.1
p_IDU_stop = 1 # probability of stoping IDU
p_HIV = 0.01 # probability of getting HIV after HCV
################################Functions#########################################

#Baseline treatment time( Note: it can be used in scenarios where time between treatment and diagnosis is fixed.) 
# It shows how long does it take to be treated after diagnosis and being eligible to be treated
##################################################

#Baseline spontaneous recovery time
###################################

sr.fun = function() 
{
  logisticfunc <- function(min, max, infl, slo) {min + ((max - min)/(1 + (teta / infl)^ slo))}
  
  sr.t <- rep(0, cohortSize)
  for (i in 1:cohortSize)
  {
    teta = runif( 1, 0, 1)  
    
    p = logisticfunc( min = 0, max = 1, infl = 0.25, slo = 2.23)
    
    w <- runif(1, 0, 1)
    ifelse( w <= p, sr.t[i] <- teta, sr.t[i] <- 999) 
  }
  sr.t
  
}#End Of Spontaneous Recovery

#####################################################################################
#BACKGROUND MORTALITY
# will be used in the hazard_fun file to calculate the harzard function for mortality. 
#To match age category with death rate (cuts is the value of the start of the intervall eg for eg 20-24 => cuts=20)
#########################################################

pw.eval.ext = function(Cuts, x, func.values, tol = 0.0001) 
{
  
  lengthCuts = length(Cuts)
  unlist(lapply(x,
                function(xp) {
                  
                  if (xp <  max(Cuts))
                    func.values[ which((Cuts[c(1:(lengthCuts - 1))] - rep(xp, lengthCuts - 1))*
                                         (Cuts[c(2:lengthCuts)] - tol  - rep(xp, lengthCuts-1)) <= 0 , arr.ind =  TRUE)[1] ]
                  
                  else func.values[lengthCuts]
                  
                }
                
                
  )) 
  
}

###################################################################################

bl <- matrix (nrow = cohortSize, ncol = bl_number)
{
  colnames(bl) <- c("Alcohol","Gender","HIV","HIV_time","MSM","IDU","IDU_time_start","IDU_time_stop", "Genotype" , "Origin", "Age", "BirthYear", "sr.time")
  ################################################################################
  ################################################################################################
  bl[,"Alcohol"] =  rep(data [case,"Alcohol"], cohortSize)
  bl[,"Genotype"] =  rep(data [case,"Genotype"], cohortSize)
  bl[,"Origin"] = rep(data [case,"Origin"], cohortSize)
  bl[,"Gender"] = rep(data[case,"Gender"], cohortSize)
  bl[,"HIV"] = rep(data[case,"HIV"], cohortSize)
  bl[,"IDU"] = rep(data[case,"IDU"], cohortSize)
  bl[,"MSM"] = rep(data[case,"MSM"], cohortSize)
  #bl[,"MSM"] = rep(1, cohortSize)
  ########################################################
  ########################################################
  age_ID = data [case,"Age"]
  bl[,"Age"] = switch (age_ID,
                       runif(cohortSize, 10, 21),
                       runif(cohortSize, 21, 31),
                       runif(cohortSize, 31, 41),
                       runif(cohortSize, 41, 51),
                       runif(cohortSize, 51, 61),
                       runif(cohortSize, 61, 71),
                       runif(cohortSize, 71, 81),
                       runif(cohortSize, 81, 91)
  )
  ####################################################
  birth_ID = data [case,"BirthYear"]
  bl[,"BirthYear"] = switch (birth_ID,
                             sample(1937: 1947, cohortSize, replace=T), 
                             sample(1948: 1957, cohortSize, replace=T),
                             sample( 1958: 1967, cohortSize, replace=T),
                             sample(1968: 1977, cohortSize, replace=T),
                             sample( 1978: 1987, cohortSize, replace=T),
                             sample(1988: 1997, cohortSize, replace=T),
                             sample(1998: 2007, cohortSize, replace=T),
                             sample(2008: 2016, cohortSize, replace=T)
  )
  ###################################################################
  bl[,"HIV_time"] <- (bl[,"HIV"])
  #######################
  for (i in 1: cohortSize)
  {
    bl[i,"HIV_time"] <- 999 
    w = runif(1, 0, 1)
    
    if (bl[i,"HIV"] == 1)
      bl[i,"HIV_time"] = 0 # time having hiv
    
    if (bl[i,"HIV"] == 0)
    {
      ID = ifelse( w < p_HIV, w, -1)
      bl[i,"HIV_time"] = ifelse( ID == -1, 999, runif(1, 0.25, maxTime)) # time having hiv
    }
  }
  ############################################################################################
  bl[,"IDU_time_start"] <- bl[,"IDU"] 
  bl[,"IDU_time_stop"] <- bl[,"IDU"]
  for (i in 1: cohortSize)
  {
    w = runif(1, 0, 1)
    if (bl[i,"IDU"] == 0)
    {
      IDU_ID = ifelse( w < p_IDU_start, w, 0)
      bl[i,"IDU_time_start"] = ifelse( IDU_ID == 0, 999, maxTime * IDU_ID + runif(1, 0, maxTime) * rbinom(1, 1, 0.5)) # time starting idu
    }
    if (bl[i,"IDU"] == 1)
      bl[i,"IDU_time_start"] = 0 # time having hiv
    teta = runif(1, 0, 1)
    bl[i,"IDU_time_stop"] = ifelse ( bl[i,"IDU_time_start"] < maxTime && teta < p_IDU_stop, bl[i,"IDU_time_start"] + 0.2 * teta * maxTime , 999)
  }
  ##########################################################################################
  bl[,"sr.time"] <- sr.fun()
}

###################################################################################
bm_M <-c(bm[,3],rep(1,70)) # to set background mortality of patients >75 years old as = background mortality of patients 75 years old 
bm_F <-c(bm[,4],rep(1,70)) # to set background mortality of patients >75 years old as = background mortality of patients 75 years old 
