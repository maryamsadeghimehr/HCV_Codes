FP_scenario <- 1
source("index_block_fs_stage.r") # it identify that in which block and stage we are in
######################################################################################

##########################################################################################
if (FP_scenario == 1)
{
  
  ###  FIBROSIS PROGRESSION VOGEL
  fp.fun = function(t, bl, history, r)
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
    ##################################### 
    res<- 0 + 0 * t
    rr<- 0
    {
      ifelse(index_block == 1 || index_block == 2 ||index_block == 3 ||index_block == 4 || index_block == 5 , rr<- 1, rr<- 0.1 ) # detectable VL/undetectable VL
      {
        ###################
        r = switch (index_fs,
                    
                    c(1,1.16,1.33),
                    c(1,1.3,2.22),
                    c(1,1.3,2.22),
                    c(1,1.16,4)
        )
        #r = c(1,1,1.33)
        ###################
        if (index_fs == 1)
        {
          if (bl["Gender"] == 1) #male
          {
            res[(t + sum(history) + bl["Age"])< 20] <- { (rr * 0.045 * r[bl["Alcohol"] + 1]) + 0 * t } 
            res[(t + sum(history)+bl["Age"])>= 20] <- { (rr*0.037 * r[bl["Alcohol"] + 1]) + 0 * t }
            res[(t + sum(history)+bl["Age"])>= 30] <-{ (rr*0.027 * r[bl["Alcohol"] + 1]) + 0 * t }
            res[(t + sum(history)+bl["Age"])>= 40] <- { (rr * 0.099 * r[bl["Alcohol"]+1]) + 0 * t }
            res[(t + sum(history)+bl["Age"])>= 50] <- { (rr * 0.121 * r[bl["Alcohol"] + 1]) + 0 * t }  
            res[(t + sum(history)+bl["Age"])>= 60] <- { (rr * 0.138 * r[bl["Alcohol"] + 1]) + 0 * t }  
            res[(t + sum(history)+bl["Age"])>= 70] <- { (rr * 0.155 * r[bl["Alcohol"] + 1]) + 0 * t }  
            res[(t + sum(history)+bl["Age"])>= 80] <- { (rr * 0.127 * r[bl["Alcohol"] + 1]) + 0 * t }  
            
          }
          if (bl["Gender"] == 0) # female
          {
            res[(t + sum(history) + bl["Age"])< 20] <- { (rr * 0.038 * r[bl["Alcohol"] + 1]) + 0 * t } 
            res[(t + sum(history)+bl["Age"])>= 20] <- { (rr*0.031 * r[bl["Alcohol"] + 1]) + 0 * t }
            res[(t + sum(history)+bl["Age"])>= 30] <-{ (rr*0.022 * r[bl["Alcohol"] + 1]) + 0 * t }
            res[(t + sum(history)+bl["Age"])>= 40] <- { (rr * 0.082 * r[bl["Alcohol"]+1]) + 0 * t }
            res[(t + sum(history)+bl["Age"])>= 50] <- { (rr * 0.102 * r[bl["Alcohol"] + 1]) + 0 * t }  
            res[(t + sum(history)+bl["Age"])>= 60] <- { (rr * 0.115 * r[bl["Alcohol"] + 1]) + 0 * t }  
            res[(t + sum(history)+bl["Age"])>= 70] <- { (rr * 0.129 * r[bl["Alcohol"] + 1]) + 0 * t }  
            res[(t + sum(history)+bl["Age"])>= 80] <- { (rr * 0.106 * r[bl["Alcohol"] + 1]) + 0 * t }
          }
        }#(index_fs == 1)
        ######################
        if (index_fs == 2)
        {
          if (bl["Gender"] == 1)
          {
            res[(t + sum(history) + bl["Age"])< 20] <- { (rr * 0.033 * r[bl["Alcohol"] + 1]) + 0 * t } 
            res[(t + sum(history)+bl["Age"])>= 20] <- { (rr*0.027 * r[bl["Alcohol"] + 1]) + 0 * t }
            res[(t + sum(history)+bl["Age"])>= 30] <-{ (rr*0.019 * r[bl["Alcohol"] + 1]) + 0 * t }
            res[(t + sum(history)+bl["Age"])>= 40] <- { (rr * 0.072 * r[bl["Alcohol"]+1]) + 0 * t }
            res[(t + sum(history)+bl["Age"])>= 50] <- { (rr * 0.088 * r[bl["Alcohol"] + 1]) + 0 * t }  
            res[(t + sum(history)+bl["Age"])>= 60] <- { (rr * 0.100 * r[bl["Alcohol"] + 1]) + 0 * t }  
            res[(t + sum(history)+bl["Age"])>= 70] <- { (rr * 0.112 * r[bl["Alcohol"] + 1]) + 0 * t }  
            res[(t + sum(history)+bl["Age"])>= 80] <- { (rr * 0.130 * r[bl["Alcohol"] + 1]) + 0 * t }    
          }
          if (bl["Gender"] == 0)
          {
            res[(t + sum(history) + bl["Age"])< 20] <- { (rr * 0.028 * r[bl["Alcohol"] + 1]) + 0 * t } 
            res[(t + sum(history)+bl["Age"])>= 20] <- { (rr*0.022 * r[bl["Alcohol"] + 1]) + 0 * t }
            res[(t + sum(history)+bl["Age"])>= 30] <-{ (rr*0.016 * r[bl["Alcohol"] + 1]) + 0 * t }
            res[(t + sum(history)+bl["Age"])>= 40] <- { (rr * 0.060 * r[bl["Alcohol"]+1]) + 0 * t }
            res[(t + sum(history)+bl["Age"])>= 50] <- { (rr * 0.074 * r[bl["Alcohol"] + 1]) + 0 * t }  
            res[(t + sum(history)+bl["Age"])>= 60] <- { (rr * 0.083 * r[bl["Alcohol"] + 1]) + 0 * t }  
            res[(t + sum(history)+bl["Age"])>= 70] <- { (rr * 0.094 * r[bl["Alcohol"] + 1]) + 0 * t }  
            res[(t + sum(history)+bl["Age"])>= 80] <- { (rr * 0.077 * r[bl["Alcohol"] + 1]) + 0 * t }  
          }
        }#(index_fs == 2)
        ######################
        if (index_fs == 3)
        {
          if (bl["Gender"] == 1)
          {
            res[(t + sum(history) + bl["Age"])< 20] <- { (rr * 0.047 * r[bl["Alcohol"] + 1]) + 0 * t } 
            res[(t + sum(history)+bl["Age"])>= 20] <- { (rr*0.038 * r[bl["Alcohol"] + 1]) + 0 * t }
            res[(t + sum(history)+bl["Age"])>= 30] <-{ (rr*0.028 * r[bl["Alcohol"] + 1]) + 0 * t }
            res[(t + sum(history)+bl["Age"])>= 40] <- { (rr * 0.102 * r[bl["Alcohol"]+1]) + 0 * t }
            res[(t + sum(history)+bl["Age"])>= 50] <- { (rr * 0.124 * r[bl["Alcohol"] + 1]) + 0 * t }  
            res[(t + sum(history)+bl["Age"])>= 60] <- { (rr * 0.141 * r[bl["Alcohol"] + 1]) + 0 * t }  
            res[(t + sum(history)+bl["Age"])>= 70] <- { (rr * 0.159 * r[bl["Alcohol"] + 1]) + 0 * t }  
            res[(t + sum(history)+bl["Age"])>= 80] <- { (rr * 0.130 * r[bl["Alcohol"] + 1]) + 0 * t }  
          }
          if (bl["Gender"] == 0)
          {
            res[(t + sum(history) + bl["Age"])< 20] <- { (rr * 0.039 * r[bl["Alcohol"] + 1]) + 0 * t } 
            res[(t + sum(history)+bl["Age"])>= 20] <- { (rr*0.031 * r[bl["Alcohol"] + 1]) + 0 * t }
            res[(t + sum(history)+bl["Age"])>= 30] <-{ (rr*0.023 * r[bl["Alcohol"] + 1]) + 0 * t }
            res[(t + sum(history)+bl["Age"])>= 40] <- { (rr * 0.085 * r[bl["Alcohol"]+1]) + 0 * t }
            res[(t + sum(history)+bl["Age"])>= 50] <- { (rr * 0.104 * r[bl["Alcohol"] + 1]) + 0 * t }  
            res[(t + sum(history)+bl["Age"])>= 60] <- { (rr * 0.118 * r[bl["Alcohol"] + 1]) + 0 * t }  
            res[(t + sum(history)+bl["Age"])>= 70] <- { (rr * 0.132 * r[bl["Alcohol"] + 1]) + 0 * t }  
            res[(t + sum(history)+bl["Age"])>= 80] <- { (rr * 0.109 * r[bl["Alcohol"] + 1]) + 0 * t }   
          }
        }#(index_fs == 3)
        ######################
        if (index_fs == 4)
        {
          if (bl["Gender"] == 1)
          {
            res[(t + sum(history) + bl["Age"])< 20] <- { (rr * 0.006 * r[bl["Alcohol"] + 1]) + 0 * t } 
            res[(t + sum(history)+bl["Age"])>= 20] <- { (rr*0.018 * r[bl["Alcohol"] + 1]) + 0 * t }
            res[(t + sum(history)+bl["Age"])>= 30] <-{ (rr*0.040 * r[bl["Alcohol"] + 1]) + 0 * t }
            res[(t + sum(history)+bl["Age"])>= 40] <- { (rr * 0.063 * r[bl["Alcohol"]+1]) + 0 * t }
            res[(t + sum(history)+bl["Age"])>= 50] <- { (rr * 0.034 * r[bl["Alcohol"] + 1]) + 0 * t }  
            res[(t + sum(history)+bl["Age"])>= 60] <- { (rr * 0.070 * r[bl["Alcohol"] + 1]) + 0 * t }  
            res[(t + sum(history)+bl["Age"])>= 70] <- { (rr * 0.136 * r[bl["Alcohol"] + 1]) + 0 * t }  
            res[(t + sum(history)+bl["Age"])>= 80] <- { (rr * 0.136 * r[bl["Alcohol"] + 1]) + 0 * t }  
          }
          if (bl["Gender"] == 0)
          {
            res[(t + sum(history) + bl["Age"])< 20] <- { (rr * 0.004 * r[bl["Alcohol"] + 1]) + 0 * t } 
            res[(t + sum(history)+bl["Age"])>= 20] <- { (rr*0.015 * r[bl["Alcohol"] + 1]) + 0 * t }
            res[(t + sum(history)+bl["Age"])>= 30] <-{ (rr*0.033 * r[bl["Alcohol"] + 1]) + 0 * t }
            res[(t + sum(history)+bl["Age"])>= 40] <- { (rr * 0.053 * r[bl["Alcohol"]+1]) + 0 * t }
            res[(t + sum(history)+bl["Age"])>= 50] <- { (rr * 0.028 * r[bl["Alcohol"] + 1]) + 0 * t }  
            res[(t + sum(history)+bl["Age"])>= 60] <- { (rr * 0.059 * r[bl["Alcohol"] + 1]) + 0 * t }  
            res[(t + sum(history)+bl["Age"])>= 70] <- { (rr * 0.113 * r[bl["Alcohol"] + 1]) + 0 * t }  
            res[(t + sum(history)+bl["Age"])>= 80] <- { (rr * 0.113 * r[bl["Alcohol"] + 1]) + 0 * t }  
          }
        }#(index_fs == 4)
        #######################
        
      }
    }
    
    return(res)
  } #  FIBROSIS PROGRESSION VOGEL
  
}# end if (FP_scenario == 1)
#######################################################################################
##########################################################################################
if (FP_scenario == 2)
{
  fp.fun = function(t, lr, bl, history){
    
    # identify the state we are in:
    m=gems:::auxcounter(ns) # matrix of all transitions (without the terminal states)
    if(all(history==0)) {index_trans=1} else {index_trans=max(which(history>0))}
    index_state=which(m==index_trans, arr.ind=TRUE)[2] # to get the colum number(arrival state) that corresponds to the last transition
    
    #return in which block is a certain state:
    ref0 = cbind(as.vector(apply(matrix(1:blocks_nb,blocks_nb,1), 1, function(x) rep(x,blocks_size))), 1:(ns-3))
    index_block.fun = function(state, ref=ref0){
      ref[which(ref[,2] == state, arr.ind = TRUE),1]
    } 
    index_block=index_block.fun(state = index_state, ref = ref0)
    
    
    #return in which fibrosis state is a certain state (which level of Fibrosis prog?)
    ref0=cbind(rep(1:blocks_size, blocks_nb), 1:(ns-3))
    index_fs.fun = function(state, ref=ref0){
      ref[which(ref[,2] == state, arr.ind = TRUE),1]
    }
    index_fs=index_fs.fun(state = index_state, ref = ref0)  
    #   res<-0+0*t
    rr<-1
    
    ifelse(index_block== 3 || index_block== 4,  rr<-0.1, rr<-1 ) #undetectable HCV (under treatment or SVR/SC) / detectable HCV VL
    
    
    res <- rep(0, length(t))  
    
    {
      res[bl["Age"]<21] <-  { 0.091*rr*exp(lr) +0*t } # we fix r=1 in params 
      res[bl["Age"]>=21] <- { 0.105*rr*exp(lr) +0*t }
      res[bl["Age"]>=31] <- { 0.138*rr*exp(lr) +0*t }
      res[bl["Age"]>=41] <- { 0.2 *rr*exp(lr) +0*t }
      res[bl["Age"]>=51] <- { 0.333*rr*exp(lr) +0*t } 
    } 
    
    return(res)
    
  }
  
}
