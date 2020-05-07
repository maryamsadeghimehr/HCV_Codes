source("index_block_fs_stage.R")
#scenario_diag = 1
scenario_diag = diag_ind
#########################################################
################################################################
if (scenario_diag == 1)
{
  ################################################################
  diag.fun = function(bl,history,t)
  {
    ## Identify in which stage and block the patient is located
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
    res <- 0 * t
    c <- 0 * t
    for (j in 1:length(t))
      ###############################################################
    if (t[j] + sum(history) + bl["Age"] + bl["BirthYear"] < 2018) #calendar time dependency
    {
      if (bl["MSM"] == 0)
        c_HIV = ifelse((bl["HIV_time"] <= t[j] + sum(history)), 0.2, 0) # for more than one year
      else
        c_HIV = 0 # for only HIV and not MSM
      
      c_IDU = ifelse((bl["IDU_time_start"] < t[j] + sum(history) && bl["IDU_time_stop"] > t[j] + sum(history)) , 0.5, 0)
      
      c_MSM_HIV = ifelse((bl["MSM"] == 1 && (bl["HIV_time"] + 1 < t[j] + sum(history))) , 1 , 0)
      
      c_MSM = ifelse((bl["MSM"] == 1 && (bl["HIV_time"] + 1 > t[j] + sum(history))) , 0.5 , 0)
      
      c_birth =  ifelse((bl["BirthYear"] >= 1951  && (bl["BirthYear"] <= 1985)) , 0 , 0)
      
      c_origin = ifelse((is.element(bl["origin"],c(1, 2, 3, 4))) , 0 , 0)
      
      c_symp  = switch(index_fs, 0, 0, 0, 1, 2, 5, 5, 5)
      c_ExIDU = ifelse((bl["IDU_time_stop"] < t[j] + sum(history)) , 0, 0) # rate that is applied to patients from 2018 after stopping IDU
      
      w <- runif(1, 0, 1)
      c_background = ifelse( w <= 1, 0.01, 0)
      c[j] = c_HIV + c_IDU + c_symp + c_birth + c_origin + c_MSM_HIV + c_MSM + c_background + c_ExIDU
    }
    else if (t[j] + sum(history) + bl["Age"] + bl["BirthYear"] >= 2018)
    {
      if (bl["MSM"] == 0)
        c_HIV = ifelse((bl["HIV_time"] <= t[j] + sum(history)), 0.2, 0) # for more than one year
      else
        c_HIV = 0 # for only HIV and not MSM
      
      c_IDU = ifelse((bl["IDU_time_start"] < t[j] + sum(history) && bl["IDU_time_stop"] > t[j] + sum(history)) , 0.5, 0)
      
      c_MSM_HIV = ifelse((bl["MSM"] == 1 && (bl["HIV_time"] + 1 < t[j] + sum(history))) , 1 , 0)
      
      c_MSM = ifelse((bl["MSM"] == 1 && (bl["HIV_time"] + 1 > t[j] + sum(history))) , 0.5 , 0)
      
      c_birth =  ifelse((bl["BirthYear"] >= 1951  && (bl["BirthYear"] <= 1985)) , 0 , 0)
      
      c_origin = ifelse((is.element(bl["origin"],c(1, 2, 3, 4))) , 0 , 0)
      
      c_symp  = switch(index_fs, 0, 0, 0, 1, 2, 5, 5, 5)
      c_ExIDU = ifelse((bl["IDU_time_stop"] < t[j] + sum(history)) , 0, 0) # rate that is applied to patients from 2018 after stopping IDU
      
      w <- runif(1, 0, 1)
      c_background = ifelse( w <= 1, 0.01, 0)
      c[j] = c_HIV + c_IDU + c_symp + c_birth + c_origin + c_MSM_HIV + c_MSM + c_background + c_ExIDU
    }
    
    
    ###########################################################################################  
    #######################################################################################
    p_rna = 1
    p_atb <- (t + sum(history) > 0.5) + (t + sum(history) <= 0.5) * 0.975 * (t + sum(history)) ^ 0.273   # we also can fit a line from 0 to 1 for the first 6 month (y = 0.494 * x + 0.45746) look at the graphs with x = c(0,0.1,0.2, 0.3,0.4,0.5,0.6,0.7, 0.8, 0.9, 1)
    # y = 0.975*(x^0.273); plot(x,y, type = "line"); lines(x,0.45*x + 0.49);
    alpha = 0.01 # the proportion of the individuals who will be tested again
    beta = 0.82 # is the proportion of the individuals who have antibody positive test 
    S.t <- 0
    ############################
    S0 <- exp(-(c * t))
    S.t = S0 + (1 - S0) * (1 - p_rna * p_atb * beta) - alpha * (1 - S0) * (1 - p_rna * p_atb * beta) * p_rna * p_atb * beta
    f.t = -c * S0 + c * S0 * (1 - p_rna * p_atb * beta) - alpha * c * S0 * (1 - p_rna * p_atb * beta) * p_rna * p_atb * beta # using the derivative of S0
    
    #####################################################################
    
    h.t = -f.t / S.t
    h.t = h.t * 0 + c # to see if a constant hazard function works better
    res <- h.t # h.t is a vector with the length of t
    for (j in 1:length(t))
      if (t[j] + sum(history) + bl["Age"] + bl["BirthYear"] < 1971) # to capture the fact that HCV has been intruduced since 1971 
        res[j] <- 0 * t[j]
    return(res)
  }
  
  
}



if (scenario_diag == 2) # for screening birth cohort
{
  ################################################################
  diag.fun = function(bl,history,t)
  {
    ## Identify in which stage and block the patient is located
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
    res <- 0 * t
    c <- 0 * t
    for (j in 1:length(t))
      ###############################################################
    if (t[j] + sum(history) + bl["Age"] + bl["BirthYear"] < 2018) #calendar time dependency
    {
      if (bl["MSM"] == 0)
        c_HIV = ifelse((bl["HIV_time"] <= t[j] + sum(history)), 0.2, 0) # for more than one year
      else
        c_HIV = 0 # for only HIV and not MSM
      
      c_IDU = ifelse((bl["IDU_time_start"] < t[j] + sum(history) && bl["IDU_time_stop"] > t[j] + sum(history)) , 0.5, 0)
      
      c_MSM_HIV = ifelse((bl["MSM"] == 1 && (bl["HIV_time"] + 1 < t[j] + sum(history))) , 1 , 0)
      
      c_MSM = ifelse((bl["MSM"] == 1 && (bl["HIV_time"] + 1 > t[j] + sum(history))) , 0.5 , 0)
      
      c_birth =  ifelse((bl["BirthYear"] >= 1951  && (bl["BirthYear"] <= 1985)) , 0 , 0)
      
      c_origin = ifelse((is.element(bl["origin"],c(1, 2, 3, 4))) , 0 , 0)
      
      c_symp  = switch(index_fs, 0, 0, 0, 1, 2, 5, 5, 5)
      c_ExIDU = ifelse((bl["IDU_time_stop"] < t[j] + sum(history)) , 0, 0) # rate that is applied to patients from 2018 after stopping IDU
      
      w <- runif(1, 0, 1)
      c_background = ifelse( w <= 1, 0.01, 0)
      c[j] = c_HIV + c_IDU + c_symp + c_birth + c_origin + c_MSM_HIV + c_MSM + c_background + c_ExIDU
    }
    else if (t[j] + sum(history) + bl["Age"] + bl["BirthYear"] >= 2018)
    {
      if (bl["MSM"] == 0)
        c_HIV = ifelse((bl["HIV_time"] <= t[j] + sum(history)), 0.2, 0) # for more than one year
      else
        c_HIV = 0 # for only HIV and not MSM
      
      c_IDU = ifelse((bl["IDU_time_start"] < t[j] + sum(history) && bl["IDU_time_stop"] > t[j] + sum(history)) , 0.5, 0)
      
      c_MSM_HIV = ifelse((bl["MSM"] == 1 && (bl["HIV_time"] + 1 < t[j] + sum(history))) , 1 , 0)
      
      c_MSM = ifelse((bl["MSM"] == 1 && (bl["HIV_time"] + 1 > t[j] + sum(history))) , 0.5 , 0)
      
      c_birth =  ifelse((bl["BirthYear"] >= 1951  && (bl["BirthYear"] <= 1985)) , 0.5 , 0)
      
      c_origin = ifelse((is.element(bl["origin"],c(1, 2, 3, 4))) , 0 , 0)
      
      c_symp  = switch(index_fs, 0, 0, 0, 1, 2, 5, 5, 5)
      c_ExIDU = ifelse((bl["IDU_time_stop"] < t[j] + sum(history)) , 0, 0) # rate that is applied to patients from 2018 after stopping IDU
      
      w <- runif(1, 0, 1)
      c_background = ifelse( w <= 1, 0.01, 0)
      c[j] = c_HIV + c_IDU + c_symp + c_birth + c_origin + c_MSM_HIV + c_MSM + c_background + c_ExIDU
    }
    
    
    ###########################################################################################  
    #######################################################################################
    p_rna = 1
    p_atb <- (t + sum(history) > 0.5) + (t + sum(history) <= 0.5) * 0.975 * (t + sum(history)) ^ 0.273   # we also can fit a line from 0 to 1 for the first 6 month (y = 0.494 * x + 0.45746) look at the graphs with x = c(0,0.1,0.2, 0.3,0.4,0.5,0.6,0.7, 0.8, 0.9, 1)
    # y = 0.975*(x^0.273); plot(x,y, type = "line"); lines(x,0.45*x + 0.49);
    alpha = 0.01 # the proportion of the individuals who will be tested again
    beta = 0.82 # is the proportion of the individuals who have antibody positive test 
    S.t <- 0
    ############################
    S0 <- exp(-(c * t))
    S.t = S0 + (1 - S0) * (1 - p_rna * p_atb * beta) - alpha * (1 - S0) * (1 - p_rna * p_atb * beta) * p_rna * p_atb * beta
    f.t = -c * S0 + c * S0 * (1 - p_rna * p_atb * beta) - alpha * c * S0 * (1 - p_rna * p_atb * beta) * p_rna * p_atb * beta # using the derivative of S0
    
    #####################################################################
    
    h.t = -f.t / S.t
    res <- h.t # h.t is a vector with the length of t
    for (j in 1:length(t))
      if (t[j] + sum(history) + bl["Age"] + bl["BirthYear"] < 1971) # to capture the fact that HCV has been intruduced since 1971 
        res[j] <- 0 * t[j]
    return(res)
  }
  
  
}




#######################################################################

################################################################
if (scenario_diag == 3) #intensive IDU
{
  ################################################################
  diag.fun = function(bl,history,t)
  {
    ## Identify in which stage and block the patient is located
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
    res <- 0 * t
    c <- 0 * t
    for (j in 1:length(t))
      ###############################################################
    if (t[j] + sum(history) + bl["Age"] + bl["BirthYear"] < 2018) #calendar time dependency
    {
      if (bl["MSM"] == 0)
        c_HIV = ifelse((bl["HIV_time"] <= t[j] + sum(history)), 0.2, 0) # for more than one year
      else
        c_HIV = 0 # for only HIV and not MSM
      
      c_IDU = ifelse((bl["IDU_time_start"] < t[j] + sum(history) && bl["IDU_time_stop"] > t[j] + sum(history)) , 0.5, 0)
      
      c_MSM_HIV = ifelse((bl["MSM"] == 1 && (bl["HIV_time"] + 1 < t[j] + sum(history))) , 1 , 0)
      
      c_MSM = ifelse((bl["MSM"] == 1 && (bl["HIV_time"] + 1 > t[j] + sum(history))) , 0.5 , 0)
      
      c_birth =  ifelse((bl["BirthYear"] >= 1951  && (bl["BirthYear"] <= 1985)) , 0 , 0)
      
      c_origin = ifelse((is.element(bl["origin"],c(1, 2, 3, 4))) , 0 , 0)
      
      c_symp  = switch(index_fs, 0, 0, 0, 1, 2, 5, 5, 5)
      c_ExIDU = ifelse((bl["IDU_time_stop"] < t[j] + sum(history)) , 0, 0) # rate that is applied to patients from 2018 after stopping IDU
      
      w <- runif(1, 0, 1)
      c_background = ifelse( w <= 1, 0.01, 0)
      c[j] = c_HIV + c_IDU + c_symp + c_birth + c_origin + c_MSM_HIV + c_MSM + c_background + c_ExIDU
    }
    else if (t[j] + sum(history) + bl["Age"] + bl["BirthYear"] >= 2018)
    {
      if (bl["MSM"] == 0)
        c_HIV = ifelse((bl["HIV_time"] <= t[j] + sum(history)), 0.2, 0) # for more than one year
      else
        c_HIV = 0 # for only HIV and not MSM
      
      c_IDU = ifelse((bl["IDU_time_start"] < t[j] + sum(history) && bl["IDU_time_stop"] > t[j] + sum(history)) , 1, 0)
      
      c_MSM_HIV = ifelse((bl["MSM"] == 1 && (bl["HIV_time"] + 1 < t[j] + sum(history))) , 1 , 0)
      
      c_MSM = ifelse((bl["MSM"] == 1 && (bl["HIV_time"] + 1 > t[j] + sum(history))) , 0.5 , 0)
      
      c_birth =  ifelse((bl["BirthYear"] >= 1951  && (bl["BirthYear"] <= 1985)) , 0 , 0)
      
      c_origin = ifelse((is.element(bl["origin"],c(1, 2, 3, 4))) , 0 , 0)
      
      c_symp  = switch(index_fs, 0, 0, 0, 1, 2, 5, 5, 5)
      c_ExIDU = ifelse((bl["IDU_time_stop"] < t[j] + sum(history)) , 0, 0) # rate that is applied to patients from 2018 after stopping IDU
      
      w <- runif(1, 0, 1)
      c_background = ifelse( w <= 1, 0.01, 0)
      c[j] = c_HIV + c_IDU + c_symp + c_birth + c_origin + c_MSM_HIV + c_MSM + c_background + c_ExIDU
    }
    
    
    ###########################################################################################  
    #######################################################################################
    p_rna = 1
    p_atb <- (t + sum(history) > 0.5) + (t + sum(history) <= 0.5) * 0.975 * (t + sum(history)) ^ 0.273  # we also can fit a line from 0 to 1 for the first 6 month (y = 0.494 * x + 0.45746) look at the graphs with x = c(0,0.1,0.2, 0.3,0.4,0.5,0.6,0.7, 0.8, 0.9, 1)
    # y = 0.975*(x^0.273); plot(x,y, type = "line"); lines(x,0.45*x + 0.49);
    alpha = 0.01 # the proportion of the individuals who will be tested again
    beta = 0.82 # is the proportion of the individuals who have antibody positive test 
    S.t <- 0
    ############################
    S0 <- exp(-(c * t))
    S.t = S0 + (1 - S0) * (1 - p_rna * p_atb * beta) - alpha * (1 - S0) * (1 - p_rna * p_atb * beta) * p_rna * p_atb * beta
    f.t = -c * S0 + c * S0 * (1 - p_rna * p_atb * beta) - alpha * c * S0 * (1 - p_rna * p_atb * beta) * p_rna * p_atb * beta # using the derivative of S0
    
    #####################################################################
    
    h.t = -f.t / S.t
    res <- h.t # h.t is a vector with the length of t
    for (j in 1:length(t))
      if (t[j] + sum(history) + bl["Age"] + bl["BirthYear"] < 1971) # to capture the fact that HCV has been intruduced since 1971 
        res[j] <- 0 * t[j]
    return(res)
  }
  
  
}

##################################
# if (scenario_diag == 4) # ex_IDU
# {
#   ################################################################
#   diag.fun = function(bl,history,t)
#   {
#     ## Identify in which stage and block the patient is located
#     statesNum = 52
#     blocks_nb = 6
#     blocks_size = 5
#     # identify the state we are in:
#     m = gems:::auxcounter(statesNum) # matrix of all transitions (without the terminal states)
#     if (all(history == 0)) {index_trans = 1} else {index_trans = max(which(history > 0))}
#     index_state = which(m == index_trans, arr.ind = TRUE)[2] # to get the colum number(arrival state) that corresponds to the last transition
#     ###################################################################
#     #return in which block is a certain state:
#     ref0 = cbind(as.vector(apply(matrix(1:blocks_nb,blocks_nb,1), 1, function(x) rep(x,blocks_size))), 1:(statesNum - 22))
#     ref1 = cbind(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6), c(31,37,43,32,38,44,33,39,45,34,40,46,35,41,47,36,42,48))
#     ref2 = rbind(ref0, ref1)
#     ########################################################################
#     #return in which fibrosis state is a certain state (which level of Fibrosis prog?)
#     ref0 = cbind(rep(1:blocks_size, blocks_nb), 1:(ns - 22))
#     ref1 = cbind(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6), c(31,37,43,32,38,44,33,39,45,34,40,46,35,41,47,36,42,48))
#     ref2 = rbind(ref0,ref1)
#     ############################
#     index_fs.fun = function(state, ref = ref2)
#     {
#       ref[which(ref[,2] == state, arr.ind = TRUE),1]
#     }
#     index_fs = index_fs.fun(state = index_state, ref = ref2)
#     state = index_state
# 
#     if (is.element(state, c(31:36)) )
#       index_fs = 6 # DC
# 
#     if (is.element(state, c(37:42)) )
#       index_fs = 7 # HCC
#     if (is.element(state, c(43:48)) )
#       index_fs = 8 #LT
#     res <- 0 * t
#     c <- 0 * t
#     for (j in 1:length(t))
#       ###############################################################
#     if (t[j] + sum(history) + bl["Age"] + bl["BirthYear"] < 2018) #calendar time dependency
#     {
#       if (bl["MSM"] == 0)
#         c_HIV = ifelse((bl["HIV_time"] <= t[j] + sum(history)), 0.2, 0) # for more than one year
#       else
#         c_HIV = 0 # for only HIV and not MSM
# 
#       c_IDU = ifelse((bl["IDU_time_start"] < t[j] + sum(history) && bl["IDU_time_stop"] > t[j] + sum(history)) , 0.5, 0)
# 
#       c_MSM_HIV = ifelse((bl["MSM"] == 1 && (bl["HIV_time"] + 1 < t[j] + sum(history))) , 1 , 0)
# 
#       c_MSM = ifelse((bl["MSM"] == 1 && (bl["HIV_time"] + 1 > t[j] + sum(history))) , 0.5 , 0)
# 
#       c_birth =  ifelse((bl["BirthYear"] >= 1951  && (bl["BirthYear"] <= 1985)) , 0 , 0)
# 
#       c_origin = ifelse((is.element(bl["origin"],c(1, 2, 3, 4))) , 0 , 0)
# 
#       c_symp  = switch(index_fs, 0, 0, 0, 1, 2, 5, 5, 5)
#       c_ExIDU = ifelse((bl["IDU_time_stop"] < t[j] + sum(history)) , 0, 0) # rate that is applied to patients from 2018 after stopping IDU
# 
#       w <- runif(1, 0, 1)
#       c_background = ifelse( w <= 1, 0.01, 0)
#       c[j] = c_HIV + c_IDU + c_symp + c_birth + c_origin + c_MSM_HIV + c_MSM + c_background + c_ExIDU
#     }
#     else if (t[j] + sum(history) + bl["Age"] + bl["BirthYear"] >= 2018)
#     {
#       if (bl["MSM"] == 0)
#         c_HIV = ifelse((bl["HIV_time"] <= t[j] + sum(history)), 0.2, 0) # for more than one year
#       else
#         c_HIV = 0 # for only HIV and not MSM
# 
#       c_IDU = ifelse((bl["IDU_time_start"] < t[j] + sum(history) && bl["IDU_time_stop"] > t[j] + sum(history)) , 0.5, 0)
# 
#       c_MSM_HIV = ifelse((bl["MSM"] == 1 && (bl["HIV_time"] + 1 < t[j] + sum(history))) , 1 , 0)
# 
#       c_MSM = ifelse((bl["MSM"] == 1 && (bl["HIV_time"] + 1 > t[j] + sum(history))) , 0.5 , 0)
# 
#       c_birth =  ifelse((bl["BirthYear"] >= 1951  && (bl["BirthYear"] <= 1985)) , 0 , 0)
# 
#       c_origin = ifelse((is.element(bl["origin"],c(1, 2, 3, 4))) , 0 , 0)
# 
#       c_symp  = switch(index_fs, 0, 0, 0, 1, 2, 5, 5, 5)
#       c_ExIDU = ifelse((bl["IDU_time_stop"] < t[j] + sum(history)) , 0.5, 0) # rate that is applied to patients from 2018 after stopping IDU
# 
#             w <- runif(1, 0, 1)
#       c_background = ifelse( w <= 1, 0.01, 0)
#       c[j] = c_HIV + c_IDU + c_symp + c_birth + c_origin + c_MSM_HIV + c_MSM + c_background + c_ExIDU
#     }
# 
#     ###########################################################################################
#     #######################################################################################
#     p_rna = 1
#     p_atb <- (t + sum(history) > 0.5) + (t + sum(history) <= 0.5) * 0.975 * (t + sum(history)) ^ 0.273  # we also can fit a line from 0 to 1 for the first 6 month (y = 0.494 * x + 0.45746) look at the graphs with x = c(0,0.1,0.2, 0.3,0.4,0.5,0.6,0.7, 0.8, 0.9, 1)
#     # y = 0.975*(x^0.273); plot(x,y, type = "line"); lines(x,0.45*x + 0.49);
#     alpha = 0.01 # the proportion of the individuals who will be tested again
#     beta = 0.82 # is the proportion of the individuals who have antibody positive test
#     S.t <- 0
#     ############################
#     S0 <- exp(-(c * t))
#     S.t = S0 + (1 - S0) * (1 - p_rna * p_atb * beta) - alpha * (1 - S0) * (1 - p_rna * p_atb * beta) * p_rna * p_atb * beta
#     f.t = -c * S0 + c * S0 * (1 - p_rna * p_atb * beta) - alpha * c * S0 * (1 - p_rna * p_atb * beta) * p_rna * p_atb * beta # using the derivative of S0
# 
#     #####################################################################
# 
#     h.t = -f.t / S.t
#     res <- h.t # h.t is a vector with the length of t
#     for (j in 1:length(t))
#       if (t[j] + sum(history) + bl["Age"] + bl["BirthYear"] < 1971) # to capture the fact that HCV has been intruduced since 1971
#         res[j] <- 0 * t[j]
#     return(res)
#   }
# 
# 
# }
# 
# 
# 
# 
# 
# 

if (scenario_diag == 4) # ex_IDU
{
  ################################################################
  diag.fun = function(bl,history,t)
  {
    ## Identify in which stage and block the patient is located
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
    res <- 0 * t
    c <- 0 * t
    ###############################################################
    
    c_HIV = ifelse((bl["HIV_time"] <= t + sum(history)&bl["MSM"] == 0), 0.2, 0) # for more than one year
    c_IDU = ifelse((bl["IDU_time_start"] < t + sum(history) && bl["IDU_time_stop"] > t + sum(history)) , 0.5, 0)
    c_MSM_HIV = ifelse((bl["MSM"] == 1 && (bl["HIV_time"] + 1 < t+ sum(history))) , 1 , 0)
    
    c_MSM = ifelse((bl["MSM"] == 1 && (bl["HIV_time"] + 1 > t + sum(history))) , 0.5 , 0)
    
    c_birth =  ifelse((bl["BirthYear"] >= 1951  && (bl["BirthYear"] <= 1985)) , 0 , 0)
    
    c_origin = ifelse((is.element(bl["origin"],c(1, 2, 3, 4))) , 0 , 0)
    
    c_symp  = switch(index_fs, 0, 0, 0, 1, 2, 5, 5, 5)
    c_ExIDU = ifelse((bl["IDU_time_stop"] < t + sum(history)) &(t + sum(history) + bl["Age"] + bl["BirthYear"] >= 2018), 0, 0) # rate that is applied to patients from 2018 after stopping IDU
    
    w <- runif(1, 0, 1)
    c_background = ifelse( w <= 1, 0.01, 0)
    c = c_HIV + c_IDU + c_symp + c_birth + c_origin + c_MSM_HIV + c_MSM + c_background + c_ExIDU
    
    ###########################################################################################
    #######################################################################################
    p_rna = 1
    p_atb <- (t + sum(history) > 0.5) + (t + sum(history) <= 0.5) * 0.975 * (t + sum(history)) ^ 0.273  # we also can fit a line from 0 to 1 for the first 6 month (y = 0.494 * x + 0.45746) look at the graphs with x = c(0,0.1,0.2, 0.3,0.4,0.5,0.6,0.7, 0.8, 0.9, 1)
    # y = 0.975*(x^0.273); plot(x,y, type = "line"); lines(x,0.45*x + 0.49);
    alpha = 0.01 # the proportion of the individuals who will be tested again
    beta = 0.82 # is the proportion of the individuals who have antibody positive test
    S.t <- 0
    ############################
    S0 <- exp(-(c * t))
    S.t = S0 + (1 - S0) * (1 - p_rna * p_atb * beta) - alpha * (1 - S0) * (1 - p_rna * p_atb * beta) * p_rna * p_atb * beta
    f.t = -c * S0 + c * S0 * (1 - p_rna * p_atb * beta) - alpha * c * S0 * (1 - p_rna * p_atb * beta) * p_rna * p_atb * beta # using the derivative of S0
    #####################################################################
    zz= ifelse(bl["IDU_time_stop"] < t + sum(history) & (t + sum(history) + bl["Age"] + bl["BirthYear"] >= 2018),0.5,0)
    h.t = -f.t / S.t
    res <- h.t+zz # h.t is a vector with the length of t
    # for (j in 1:length(t))
    #   if (t[j] + sum(history) + bl["Age"] + bl["BirthYear"] < 1971) # to capture the fact that HCV has been intruduced since 1971
    #     res[j] <- 0 * t[j]
    # return(res)
    
    res = ifelse (t + sum(history) + bl["Age"] + bl["BirthYear"] < 1971, 0 * res,res) # to capture the fact that HCV has been intruduced since 1971
    return(res)
  }
  
}

