######################################################################
### File containing all the harzard functions
### It contains fibrosis progression, spontan-
### eous recovery, diagnosis, treatment, SVR
### 
### 
#####################################################################
source("index_block_fs_stage.R") # it identify that in which block and stage we are in

#####################################################################
# Function name: fp.fun
#Description: This function calculate the progression between F0 to cirrhosis (F4)
####################################################################
source("FP_Scenarios.R")
#####################################################################
# Function name: DC.fun
#Description: This function calculate the progression from F4 to DC

dc.fun <- function(t,bl,history, lr = 0)
{
  
  statesNum = 52
  blocks_nb = 6
  blocks_size = 5
  # identify the state we are in:
  m = gems:::auxcounter(statesNum) # matrix of all transitions (without the terminal states)
  if (all(history == 0)) {index_trans = 1} else {index_trans = max(which(history > 0))}
  index_state = which(m == index_trans, arr.ind = TRUE)[2] # to get the colum number(arrival state) that corresponds to the last transition
  ###################################################################
  #return in which block is a certain state:
  ref0 = cbind(as.vector(apply(matrix(1:blocks_nb,blocks_nb,1), 1, function(x) rep(x,blocks_size))), 1:(statesNum - 22))
  ref1 = cbind(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6), c(31,37,43,32,38,44,33,39,45,34,40,46,35,41,47,36,42,48))
  ref2 = rbind(ref0, ref1)
  ########################################################################
  index_block.fun = function(state, ref=ref2)
  {
    ref[which(ref[,2] == state, arr.ind = TRUE),1]
  } 
  index_block = index_block.fun(state = index_state, ref = ref0)
  #######################################################################
  #return in which fibrosis state is a certain state (which level of Fibrosis prog?)
  ref0 = cbind(rep(1:blocks_size, blocks_nb), 1:(ns - 22))
  ref1 = cbind(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6), c(31,37,43,32,38,44,33,39,45,34,40,46,35,41,47,36,42,48))
  ref2 = rbind(ref0,ref1)
  ############################
  index_fs.fun = function(state, ref = ref2)
  {
    ref[which(ref[,2] == state, arr.ind = TRUE),1]
  }
  index_fs = index_fs.fun(state = index_state, ref = ref2)  
  state = index_state
  
  if (is.element(state, c(31:36)) )
    index_fs = 6 # DC
  
  if (is.element(state, c(37:42)) )
    index_fs = 7 # HCC
  if (is.element(state, c(43:48)) )
    index_fs = 8 #LT
  ####################################
  ifelse(index_block == 6 ,  rr <- 0.1, rr <- 1 ) # undetectable/detectable HCV VL
  r = c(1,1,1)
  ########################################################
  res <- 0 + 0 * t
  ######################################################
  if (index_fs == 5)
  {
    res[(t + sum(history) + bl["Age"]) <  30] <- {(rr * 0.0651 * r[bl["Alcohol"] + 1]) + 0 * t } 
    res[(t + sum(history) + bl["Age"]) >= 30] <- {(rr*0.0641 * r[bl["Alcohol"] + 1]) + 0 * t }
    res[(t + sum(history) + bl["Age"]) >= 40] <- {(rr*0.0648 * r[bl["Alcohol"] + 1]) + 0 * t }
    res[(t + sum(history) + bl["Age"]) >= 50] <- {(rr * 0.0649 * r[bl["Alcohol"] + 1]) + 0 * t }  
    res[(t + sum(history) + bl["Age"]) >= 60] <- {(rr * 0.0635 * r[bl["Alcohol"] + 1]) + 0 * t }  
    res[(t + sum(history) + bl["Age"]) >= 70] <- {(rr * 0.0630 * r[bl["Alcohol"] + 1]) + 0 * t }  
  }#(index_fs == 5)
  
  return(res)
}

###################################################################

# Function name: HCC.fun
#Description: This function calculate the progression from F4 to HCC

hcc.fun <- function(t,bl,history, lr)
{
  # lr = 0
  # index = stage_identifier(lr, bl, history)
  # index_block = index[1]
  # index_fs = index[2]
  # state = index[3]
  statesNum = 52
  blocks_nb = 6
  blocks_size = 5
  # identify the state we are in:
  m = gems:::auxcounter(statesNum) # matrix of all transitions (without the terminal states)
  if (all(history == 0)) {index_trans = 1} else {index_trans = max(which(history > 0))}
  index_state = which(m == index_trans, arr.ind = TRUE)[2] # to get the colum number(arrival state) that corresponds to the last transition
  ###################################################################
  #return in which block is a certain state:
  ref0 = cbind(as.vector(apply(matrix(1:blocks_nb,blocks_nb,1), 1, function(x) rep(x,blocks_size))), 1:(statesNum - 22))
  ref1 = cbind(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6), c(31,37,43,32,38,44,33,39,45,34,40,46,35,41,47,36,42,48))
  ref2 = rbind(ref0, ref1)
  ########################################################################
  index_block.fun = function(state, ref=ref2)
  {
    ref[which(ref[,2] == state, arr.ind = TRUE),1]
  } 
  index_block = index_block.fun(state = index_state, ref = ref0)
  #######################################################################
  #return in which fibrosis state is a certain state (which level of Fibrosis prog?)
  ref0 = cbind(rep(1:blocks_size, blocks_nb), 1:(ns - 22))
  ref1 = cbind(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6), c(31,37,43,32,38,44,33,39,45,34,40,46,35,41,47,36,42,48))
  ref2 = rbind(ref0,ref1)
  ############################
  index_fs.fun = function(state, ref = ref2)
  {
    ref[which(ref[,2] == state, arr.ind = TRUE),1]
  }
  index_fs = index_fs.fun(state = index_state, ref = ref2)  
  state = index_state
  
  if (is.element(state, c(31:36)) )
    index_fs = 6 # DC
  
  if (is.element(state, c(37:42)) )
    index_fs = 7 # HCC
  if (is.element(state, c(43:48)) )
    index_fs = 8 #LT
  
  #####################################################
  ifelse(index_block == 6 ,  rr <- 0.1, rr <- 1 ) # undetectable/detectable HCV VL
  r = c(1,1,1)
  res <- 0 + 0 * t
  #####################################################
  if (index_fs == 5)
  {
    res[(t + sum(history) + bl["Age"]) < 30] <- {(rr * 0.0079 * r[bl["Alcohol"] + 1]) + 0 * t } 
    res[(t + sum(history) + bl["Age"]) >= 30] <- {(rr*0.0130 * r[bl["Alcohol"] + 1]) + 0 * t }
    res[(t + sum(history) + bl["Age"]) >= 40] <- {(rr*0.0212 * r[bl["Alcohol"] + 1]) + 0 * t }
    res[(t + sum(history) + bl["Age"]) >= 50] <- {(rr * 0.0347 * r[bl["Alcohol"] + 1]) + 0 * t }  
    res[(t + sum(history) + bl["Age"]) >= 60] <- {(rr * 0.0565 * r[bl["Alcohol"] + 1]) + 0 * t }  
    res[(t + sum(history) + bl["Age"]) >= 70] <- {(rr * 0.0913 * r[bl["Alcohol"] + 1]) + 0 * t }  
  }#(index_fs == 5)
  return(res)
}
###############################################################

# Function name: dchcc.fun
#Description: This function calculate the progression from DC to HCC
dchcc.fun <- function(t,bl,history, lr)
{
  statesNum = 52
  blocks_nb = 6
  blocks_size = 5
  # identify the state we are in:
  m = gems:::auxcounter(statesNum) # matrix of all transitions (without the terminal states)
  if (all(history == 0)) {index_trans = 1} else {index_trans = max(which(history > 0))}
  index_state = which(m == index_trans, arr.ind = TRUE)[2] # to get the colum number(arrival state) that corresponds to the last transition
  ###################################################################
  #return in which block is a certain state:
  ref0 = cbind(as.vector(apply(matrix(1:blocks_nb,blocks_nb,1), 1, function(x) rep(x,blocks_size))), 1:(statesNum - 22))
  ref1 = cbind(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6), c(31,37,43,32,38,44,33,39,45,34,40,46,35,41,47,36,42,48))
  ref2 = rbind(ref0, ref1)
  ########################################################################
  index_block.fun = function(state, ref=ref2)
  {
    ref[which(ref[,2] == state, arr.ind = TRUE),1]
  } 
  index_block = index_block.fun(state = index_state, ref = ref0)
  #######################################################################
  #return in which fibrosis state is a certain state (which level of Fibrosis prog?)
  ref0 = cbind(rep(1:blocks_size, blocks_nb), 1:(ns - 22))
  ref1 = cbind(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6), c(31,37,43,32,38,44,33,39,45,34,40,46,35,41,47,36,42,48))
  ref2 = rbind(ref0,ref1)
  ############################
  index_fs.fun = function(state, ref = ref2)
  {
    ref[which(ref[,2] == state, arr.ind = TRUE),1]
  }
  index_fs = index_fs.fun(state = index_state, ref = ref2)  
  state = index_state
  
  if (is.element(state, c(31:36)) )
    index_fs = 6 # DC
  
  if (is.element(state, c(37:42)) )
    index_fs = 7 # HCC
  if (is.element(state, c(43:48)) )
    index_fs = 8 #LT
  #####################################################
  r = c(1,1,1)
  res <- rep(exp(lr), length(t))
  rr <- ifelse(index_block == 1 || index_block == 2 || index_block == 3 || index_block == 4 || index_block == 5 ,
               1, 0.1 ) # detectable VL/undetectable VL
  if (index_fs == 6)
  {
    res[(t + sum(history) + bl["Age"]) < 30] <- {(rr * 0.0155 * r[bl["Alcohol"] + 1]) + 0 * t } 
    res[(t + sum(history) + bl["Age"]) >= 30] <- {(rr*0.0252 * r[bl["Alcohol"] + 1]) + 0 * t }
    res[(t + sum(history) + bl["Age"]) >= 40] <- {(rr*0.0410 * r[bl["Alcohol"] + 1]) + 0 * t }
    res[(t + sum(history) + bl["Age"]) >= 50] <- {(rr * 0.0665 * r[bl["Alcohol"] + 1]) + 0 * t }  
    res[(t + sum(history) + bl["Age"]) >= 60] <- {(rr * 0.1091 * r[bl["Alcohol"] + 1]) + 0 * t }  
    res[(t + sum(history) + bl["Age"]) >= 70] <- {(rr * 0.1762 * r[bl["Alcohol"] + 1]) + 0 * t }  
  }
  return(res)
}
##############################################################

# Function name: dcLT.fun
#Description: This function calculate the progression from DC to LT
dcLT.fun <- function(t,lr, bl, history)
{
  # index = stage_identifier(lr, bl, history)
  # index_block = index[1]
  # index_fs = index[2]
  # state = index[3]
  #######################################
  statesNum = 52
  blocks_nb = 6
  blocks_size = 5
  # identify the state we are in:
  m = gems:::auxcounter(statesNum) # matrix of all transitions (without the terminal states)
  if (all(history == 0)) {index_trans = 1} else {index_trans = max(which(history > 0))}
  index_state = which(m == index_trans, arr.ind = TRUE)[2] # to get the colum number(arrival state) that corresponds to the last transition
  ###################################################################
  #return in which block is a certain state:
  ref0 = cbind(as.vector(apply(matrix(1:blocks_nb,blocks_nb,1), 1, function(x) rep(x,blocks_size))), 1:(statesNum - 22))
  ref1 = cbind(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6), c(31,37,43,32,38,44,33,39,45,34,40,46,35,41,47,36,42,48))
  ref2 = rbind(ref0, ref1)
  #return in which fibrosis state is a certain state (which level of Fibrosis prog?)
  ref0 = cbind(rep(1:blocks_size, blocks_nb), 1:(ns - 22))
  ref1 = cbind(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6), c(31,37,43,32,38,44,33,39,45,34,40,46,35,41,47,36,42,48))
  ref2 = rbind(ref0,ref1)
  ############################
  index_fs.fun = function(state, ref = ref2)
  {
    ref[which(ref[,2] == state, arr.ind = TRUE),1]
  }
  index_fs = index_fs.fun(state = index_state, ref = ref2)  
  state = index_state
  
  if (is.element(state, c(31:36)) )
    index_fs = 6 # DC
  
  if (is.element(state, c(37:42)) )
    index_fs = 7 # HCC
  if (is.element(state, c(43:48)) )
    index_fs = 8 #LT
  ##########################################
  treated_set = c(34, 35, 36)
  if (is.element(state, treated_set))
    res <- rep(exp(lr), length(t))
  else 
    res <- rep(0, length(t)) 
  return(res)
}

#############################################################

# Function name: hccLT.fun
#Description: This function calculate the progression from HCC to LT
hccLT.fun <- function(t,lr,bl,history)
{
  statesNum = 52
  blocks_nb = 6
  blocks_size = 5
  # identify the state we are in:
  m = gems:::auxcounter(statesNum) # matrix of all transitions (without the terminal states)
  if (all(history == 0)) {index_trans = 1} else {index_trans = max(which(history > 0))}
  index_state = which(m == index_trans, arr.ind = TRUE)[2] # to get the colum number(arrival state) that corresponds to the last transition
  ###################################################################
  #return in which block is a certain state:
  ref0 = cbind(as.vector(apply(matrix(1:blocks_nb,blocks_nb,1), 1, function(x) rep(x,blocks_size))), 1:(statesNum - 22))
  ref1 = cbind(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6), c(31,37,43,32,38,44,33,39,45,34,40,46,35,41,47,36,42,48))
  ref2 = rbind(ref0, ref1)
  ########################################################################
  #return in which fibrosis state is a certain state (which level of Fibrosis prog?)
  ref0 = cbind(rep(1:blocks_size, blocks_nb), 1:(ns - 22))
  ref1 = cbind(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6), c(31,37,43,32,38,44,33,39,45,34,40,46,35,41,47,36,42,48))
  ref2 = rbind(ref0,ref1)
  ############################
  index_fs.fun = function(state, ref = ref2)
  {
    ref[which(ref[,2] == state, arr.ind = TRUE),1]
  }
  index_fs = index_fs.fun(state = index_state, ref = ref2)  
  state = index_state
  
  if (is.element(state, c(31:36)) )
    index_fs = 6 # DC
  
  if (is.element(state, c(37:42)) )
    index_fs = 7 # HCC
  if (is.element(state, c(43:48)) )
    index_fs = 8 #LT
  ##########################################
  treated_set = c(40, 41, 42)
  if (is.element(state, treated_set))
    res <-  rep(exp(lr), length(t))
  else 
    res <- rep(0, length(t)) 
  return(res)
  
}

####################################### End of the vertically progression #####################################################
###############################################################################################################################



#Description: Casecade of care
###################################################

# Function name: sr.fun
# Description: This function calculate the spontaneous recovery rate 
# Note: time to spontaneous clearance is calculated before the simulation and assigned as a baseline
##################################################
source("spontaneous_recovery.R")
#####################################################################

# Function name: aToC
# Description: This function give us the fixed time from acute to chronic states
aToC.fun = function(history)
{
  chronic.t <- (6/12) - sum(history)
} # End of aToC function

###################################################################

# Diagnosis: time of diagnosis calculated as a baseline
# In order to check for the different screening strategy please run:
source("diag_scenario.R")
##########################################################
##################################################################
# Treatment uptake 
# can be used in scenarios where time between treatment and diagnosis is a baseline
#ans = "Fg"
source(sprintf("scenario_treat%s.R",ans))
################################################################
# In order to check for the different scenarios to define SVR function please run:
# source(SVR_scenario.R)
##########################################################
#SVR
svr.fun <- function(bl,history,lr = 0)
{
  statesNum = 52
  m = gems:::auxcounter(statesNum) # matrix of all transitions (without the terminal states)
  if (all(history == 0)) {index_trans = 1} else {index_trans = max(which(history > 0))}
  index_state = which(m == index_trans, arr.ind = TRUE)[2] # to get the colum number(arrival state) that corresponds to the last transition
  index_previous_state = which(m == index_trans, arr.ind = TRUE)[1] # to get the row number(departure state) that corresponds to the last transition
  history_tmp = setdiff(which(history > 0), max(which(history > 0)))
  
  w = runif(1,0,1)
  p = 0.98
  ##########################################################################
  start_treatment0 = sum(history[1:660]) + sum(history[1096:1155]) + sum(history[1207:1248]) + sum(history[1282:1305])
  ## for untreated individuals + DC + HCC + LT
  res <- ifelse(w <= p  && sum(history) - start_treatment0 < 12/52, 12/52 - (sum(history) - start_treatment0), 999)
  
  #############################################################################################
  start_treatment1 = sum(history[1:830]) + sum(history[1096:1173]) + sum(history[1207:1260]) + sum(history[1282:1311])
  ## for untreated individuals + DC + HCC + LT
  if (is.element(index_state, c(21, 22, 23, 24, 25, 35, 41, 47)))
    res <- ifelse(w <= p  && sum(history) - start_treatment1 < 12/52, 12/52 - (sum(history) - start_treatment1), 999)  
  ##########################################################################
  return(res)
  
}
######################################################################
failed_svr.fun <- function(bl,history,lr = 0)
{
  statesNum = 52
  m = gems:::auxcounter(statesNum) # matrix of all transitions (without the terminal states)
  if (all(history == 0)) {index_trans = 1} else {index_trans = max(which(history > 0))}
  
  p = 0.5   # probability of going for the second round of treatment
  w = runif(1,0,1)
  start_treatment0 = sum(history[1:660]) + sum(history[1096:1155]) + sum(history[1207:1248]) + sum(history[1282:1305])
  ## for untreated individuals + DC + HCC + LT
  res <- ifelse(w <= p  && sum(history) - start_treatment0 < 12/52 + 0.24, 12/52 - (sum(history) - start_treatment0) + 0.24, 999)
  return (res)
}

###################################################################

# Function name: mort.fun
# Description: Mortality function
mort.fun = function(t, bl, history)
{
  # lr = 0
  # # identify the state we are in:
  # index.fun = stage_identifier(lr, bl, history)
  # index_block=index.fun[1]
  # index_fs = index.fun[2]
  # index_state = index.fun[3]
  # #######################################
  # state = index.fun[3]
  #######################################
  statesNum = 52
  blocks_nb = 6
  blocks_size = 5
  # identify the state we are in:
  m = gems:::auxcounter(statesNum) # matrix of all transitions (without the terminal states)
  if (all(history == 0)) {index_trans = 1} else {index_trans = max(which(history > 0))}
  index_state = which(m == index_trans, arr.ind = TRUE)[2] # to get the colum number(arrival state) that corresponds to the last transition
  ###################################################################
  #return in which block is a certain state:
  ref0 = cbind(as.vector(apply(matrix(1:blocks_nb,blocks_nb,1), 1, function(x) rep(x,blocks_size))), 1:(statesNum - 22))
  ref1 = cbind(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6), c(31,37,43,32,38,44,33,39,45,34,40,46,35,41,47,36,42,48))
  ref2 = rbind(ref0, ref1)
  ########################################################################
  index_block.fun = function(state, ref=ref2)
  {
    ref[which(ref[,2] == state, arr.ind = TRUE),1]
  } 
  index_block = index_block.fun(state = index_state, ref = ref0)
  #######################################################################
  #return in which fibrosis state is a certain state (which level of Fibrosis prog?)
  ref0 = cbind(rep(1:blocks_size, blocks_nb), 1:(ns - 22))
  ref1 = cbind(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6), c(31,37,43,32,38,44,33,39,45,34,40,46,35,41,47,36,42,48))
  ref2 = rbind(ref0,ref1)
  ############################
  index_fs.fun = function(state, ref = ref2)
  {
    ref[which(ref[,2] == state, arr.ind = TRUE),1]
  }
  index_fs = index_fs.fun(state = index_state, ref = ref2)  
  state = index_state
  
  if (is.element(state, c(31:36)) )
    index_fs = 6 # DC
  
  if (is.element(state, c(37:42)) )
    index_fs = 7 # HCC
  if (is.element(state, c(43:48)) )
    index_fs = 8 #LT
  
  #####################################################
  gender = bl["Gender"]
  res = 0
  if (gender == 1)
    res = (pw.eval.ext(Cuts = bm[,2], x = bl["Age"] + sum(history) + t, func.values = bm_M ))
  else
    res = (pw.eval.ext(Cuts = bm[,2], x = bl["Age"] + sum(history) + t, func.values = bm_F ))
  
  if (t + sum(history) + bl["Age"] > 100)  res = rep(1000, length(t))
  
  return(res)
  
}

####################################################################
mort_IDU.fun = function(t, bl, history)
{
  # identify the state we are in:
  # index.fun = stage_identifier(lr, bl, history)
  # index_block = index.fun[1]
  # index_fs = index.fun[2]
  # index_state = index.fun[3]
  # state = index.fun[3]
  #######################################
  #######################################
  statesNum = 52
  blocks_nb = 6
  blocks_size = 5
  # identify the state we are in:
  m = gems:::auxcounter(statesNum) # matrix of all transitions (without the terminal states)
  if (all(history == 0)) {index_trans = 1} else {index_trans = max(which(history > 0))}
  index_state = which(m == index_trans, arr.ind = TRUE)[2] # to get the colum number(arrival state) that corresponds to the last transition
  ###################################################################
  #return in which block is a certain state:
  ref0 = cbind(as.vector(apply(matrix(1:blocks_nb,blocks_nb,1), 1, function(x) rep(x,blocks_size))), 1:(statesNum - 22))
  ref1 = cbind(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6), c(31,37,43,32,38,44,33,39,45,34,40,46,35,41,47,36,42,48))
  ref2 = rbind(ref0, ref1)
  ########################################################################
  index_block.fun = function(state, ref=ref2)
  {
    ref[which(ref[,2] == state, arr.ind = TRUE),1]
  } 
  index_block = index_block.fun(state = index_state, ref = ref0)
  #######################################################################
  #return in which fibrosis state is a certain state (which level of Fibrosis prog?)
  ref0 = cbind(rep(1:blocks_size, blocks_nb), 1:(ns - 22))
  ref1 = cbind(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6), c(31,37,43,32,38,44,33,39,45,34,40,46,35,41,47,36,42,48))
  ref2 = rbind(ref0,ref1)
  ############################
  index_fs.fun = function(state, ref = ref2)
  {
    ref[which(ref[,2] == state, arr.ind = TRUE),1]
  }
  index_fs = index_fs.fun(state = index_state, ref = ref2)  
  state = index_state
  
  if (is.element(state, c(31:36)) )
    index_fs = 6 # DC
  
  if (is.element(state, c(37:42)) )
    index_fs = 7 # HCC
  if (is.element(state, c(43:48)) )
    index_fs = 8 #LT
  
  #####################################################
  ####################################
  res <- rep(0, length(t))
  #######################################################
  if (bl["IDU"] == 1 && bl["IDU_time_stop"] > sum(history) + t) # if using drug and not stop it yet
  {
    beta = 0.00017 # adjustment for IDU related mortality
    gender = bl["Gender"]
    
    if (gender == 1)
      res <- rep(beta, length(t))  
    else
      res <- rep(beta, length(t)) 
  }
  #####################
  if ((bl["IDU_time_start"] < sum(history) + t  && bl["IDU_time_stop"] > sum(history) + t)) # if using drug, but stop it at some point
  {
    beta = 0.00017 # adjustment for IDU related mortality
    gender = bl["Gender"]
    if (gender == 1)
      res <- rep(beta, length(t))
    else
      res <- rep(beta, length(t))
    
  }
  #####################
  return(res)
}
#####################################################################

mort_HIV.fun = function(t, bl, history)
{
  # identify the state we are in:
  # index.fun = stage_identifier(lr, bl, history)
  # index_block = index.fun[1]
  # index_fs = index.fun[2]
  # index_state = index.fun[3]
  # state = index.fun[3]
  #######################################
  #######################################
  statesNum = 52
  blocks_nb = 6
  blocks_size = 5
  # identify the state we are in:
  m = gems:::auxcounter(statesNum) # matrix of all transitions (without the terminal states)
  if (all(history == 0)) {index_trans = 1} else {index_trans = max(which(history > 0))}
  index_state = which(m == index_trans, arr.ind = TRUE)[2] # to get the colum number(arrival state) that corresponds to the last transition
  ###################################################################
  #return in which block is a certain state:
  ref0 = cbind(as.vector(apply(matrix(1:blocks_nb,blocks_nb,1), 1, function(x) rep(x,blocks_size))), 1:(statesNum - 22))
  ref1 = cbind(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6), c(31,37,43,32,38,44,33,39,45,34,40,46,35,41,47,36,42,48))
  ref2 = rbind(ref0, ref1)
  #######################################################################
  #return in which fibrosis state is a certain state (which level of Fibrosis prog?)
  ref0 = cbind(rep(1:blocks_size, blocks_nb), 1:(ns - 22))
  ref1 = cbind(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6), c(31,37,43,32,38,44,33,39,45,34,40,46,35,41,47,36,42,48))
  ref2 = rbind(ref0,ref1)
  ############################
  index_fs.fun = function(state, ref = ref2)
  {
    ref[which(ref[,2] == state, arr.ind = TRUE),1]
  }
  index_fs = index_fs.fun(state = index_state, ref = ref2)  
  state = index_state
  
  if (is.element(state, c(31:36)) )
    index_fs = 6 # DC
  
  if (is.element(state, c(37:42)) )
    index_fs = 7 # HCC
  if (is.element(state, c(43:48)) )
    index_fs = 8 #LT
  
  #####################################################
  #####################################################
  alpha = 1/100000
  gender = bl["Gender"]
  if (gender == 1)
    res <- ifelse(bl["HIV"] == 1 | (bl["HIV_time"] <= (t + sum(history))), rep(alpha, length(t)), rep(0, length(t)) )
  else
    res <- ifelse(bl["HIV"] == 1 | (bl["HIV_time"] <= (t + sum(history))), rep(alpha, length(t)), rep(0, length(t)) )  
  return(res)
  
}

#####################################################################
mort_Liver.fun = function(t, bl, history)
{
  # identify the state we are in:
  # index.fun = stage_identifier(lr, bl, history)
  # index_block=index.fun[1]
  # index_fs = index.fun[2]
  # index_state = index.fun[3]
  # state = index.fun[3]
  #######################################
  #######################################
  statesNum = 52
  blocks_nb = 6
  blocks_size = 5
  # identify the state we are in:
  m = gems:::auxcounter(statesNum) # matrix of all transitions (without the terminal states)
  if (all(history == 0)) {index_trans = 1} else {index_trans = max(which(history > 0))}
  index_state = which(m == index_trans, arr.ind = TRUE)[2] # to get the colum number(arrival state) that corresponds to the last transition
  ###################################################################
  #return in which block is a certain state:
  ref0 = cbind(as.vector(apply(matrix(1:blocks_nb,blocks_nb,1), 1, function(x) rep(x,blocks_size))), 1:(statesNum - 22))
  ref1 = cbind(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6), c(31,37,43,32,38,44,33,39,45,34,40,46,35,41,47,36,42,48))
  ref2 = rbind(ref0, ref1)
  #######################################################################
  #return in which fibrosis state is a certain state (which level of Fibrosis prog?)
  ref0 = cbind(rep(1:blocks_size, blocks_nb), 1:(ns - 22))
  ref1 = cbind(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6), c(31,37,43,32,38,44,33,39,45,34,40,46,35,41,47,36,42,48))
  ref2 = rbind(ref0,ref1)
  ############################
  index_fs.fun = function(state, ref = ref2)
  {
    ref[which(ref[,2] == state, arr.ind = TRUE),1]
  }
  index_fs = index_fs.fun(state = index_state, ref = ref2)  
  state = index_state
  
  if (is.element(state, c(31:36)) )
    index_fs = 6 # DC
  
  if (is.element(state, c(37:42)) )
    index_fs = 7 # HCC
  if (is.element(state, c(43:48)) )
    index_fs = 8 #LT
  
  #####################################################
  ####################################
  r_cc = 0.010
  r_dc = 0.129
  r_dc1 = 0.129
  r_hcc = 0.430
  r_hcc1 = 0.430
  r_lt = 0.160
  r_lt1 = 0.057
  
  res <- rep(0, length(t))
  res[index_state == 5 | index_state == 10 | index_state == 15 | index_state == 20 | index_state == 25 | index_state == 30 ] <- r_cc + 0 * t
  res[index_state == 31 | index_state == 32 | index_state == 33 | index_state == 34 | index_state == 35 | index_state == 36 ] <- r_dc + 0 * t
  res[(index_state == 31 | index_state == 32 | index_state == 33 | index_state == 34 | index_state == 35 | index_state == 36) && t <= 1] <- r_dc1 + 0 * t
  res[index_state == 37 | index_state == 38 | index_state == 39 | index_state == 40 | index_state == 41 | index_state == 42] <- r_hcc + 0 * t
  res[(index_state == 37 | index_state == 38 | index_state == 39 | index_state == 40 | index_state == 41 | index_state == 42) && t <= 1] <- r_hcc1 + 0 * t
  res[index_state == 43 | index_state == 44 | index_state == 45 | index_state == 46 | index_state == 47 | index_state == 48] <- r_lt + 0 * t
  res[(index_state == 43 | index_state == 44 | index_state == 45 | index_state == 46 | index_state == 47 | index_state == 48) && t <= 1] <- r_lt1 + 0 * t
  
  return(res)
  
  
}


#######################################################################
