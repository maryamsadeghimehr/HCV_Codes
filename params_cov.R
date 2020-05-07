######################################################################
### This code provides all the required parameters####################
###According to sim_cohort.R function, in order to simulation we need to define the functions and matrices below.
### 
### 
######################################################################################################
#Inputs
blocks_nb = 6 ## number of blocks of states
blocks_size = 5 ## number of states per blocks 
#######################
## terminal states (endpoints) that are not part of a block
DC <- c(31:36)
HCC <- c(37:42)
LT <- c(43:48)
Death <- (49:52)
#####################################################################################################

# General set of functions for vertical progression (F0->F1, F1->F2, F2->F3, F3-<F4, F4->DC, F4-HCC, DC->HCC, DC -> LT, HCC -> LT)
func = list(fp.fun, fp.fun, fp.fun, fp.fun, dc.fun, hcc.fun, dchcc.fun,dcLT.fun, hccLT.fun)
mort.func =  list(mort.fun, mort_IDU.fun, mort_HIV.fun, mort_Liver.fun)
###########################################################
# Initial set of parameters for vertical progression
par = list(list((1), (1), (1), (1), (1)),
           list((1), (1), (1), (1), (1)),
           list((1), (1), (1), (1), (1)),
           list((1), (1), (1), (1), (1)),
           list(1),
           list(1),
           list(1),
           list(1),
           list(1)
)
par = list(list(log(1), log(1), log(1), log(1), log(1)), #F0 -> F1
           list(log(1), log(1), log(1), log(1), log(1)), # F1-> F2
           list(log(1), log(1), log(1), log(1), log(1)), # F2 -> F3
           list(log(1), log(1), log(1), log(1), log(1)), # F3 -> F4
           list(log(0.039)), #F4 -> DC
           list(log(0.014)), #F4 -> HCC
           list(log(0.014)), #DC -> HCC
           list(log(0.031)), #DC -> LT
           list(log(0.170)) #HCC -> LT
)

length_func  =  length(func)
if (length_func != length(par)) 
  
{
  stop("You need a list of arguments per function")
}
####################################################

tf = generateHazardMatrix(statesNum) #template for transitions
pm = generateParameterMatrix(tf) #template for parameters
sigma <-  generateParameterCovarianceMatrix(pm) #Covariance matrix of parameters
######################################################
## number of first state in each block
states0 =  c(1,6,11,16,21,26)
froms   =  rep(1,blocks_nb)

if (length(states0) != length(froms)) 
  stop("Each Starting State Needs a Initial Fibrosis")

block =  list()

for (i in 1: length(states0))
{
  state0   = states0[i]
  from =  froms[i]
  statedc = DC[i]
  statehcc = HCC[i]
  stateLT = LT[i]
  ############################################################
  ### setting the structure of the blocks
  #######################################################
  #transitions to set:
  #all index for the transitions
  fs = blocks_size - (from - 1) # fibrosis state 
  
  ## indexes of states between 2 blocks
  index0 = cbind((state0):(state0 + fs - 2), (state0 + 1):(state0 + fs - 1))
  
  #adding terminal states 
  index = rbind(index0, c(fs + state0 - 1,statedc), c(fs + state0 - 1,statehcc), c(statedc, statehcc), c(statedc, stateLT), c(statehcc, stateLT))
  k = 1 #dyn counter
  
  while (k <= dim(index)[1])
  {
    
    if (index[k,1] < index[k,2])
    { 
      
      tf[[index[k,1],index[k,2]]] = func[[k]]
      pm[[index[k,1],index[k,2]]] = par[k]
      k =  k + 1
      
    }
    
    else k = k + 1  
  }
  
  #########################################################
  #Crossed transitions (horizontal)// Generating blocks
  block[[i]]  = rbind(c(index0[,1], tail(index0[,2],1)), rev( (blocks_size:1)[1:(length(index0[,1]) + 1)] ))
}


######## blocks transition matrix################
block.tran.mat = matrix(0, length(states0), length(states0))
#allowed transitions between bocks:
block.tran.mat[1,2]  = 1
block.tran.mat[1,3]  = 1
block.tran.mat[1,6]  = 1
block.tran.mat[2,3]  = 1
block.tran.mat[2,6]  = 1
block.tran.mat[3,6]  = 1
block.tran.mat[3,4]  = 1
block.tran.mat[4,5]  = 1
block.tran.mat[4,6]  = 1
block.tran.mat[5,6]  = 1
##list of lists of funtion per transition###########################
sr.func.list = list(sr.fun,sr.fun, sr.fun, sr.fun, sr.fun, sr.fun)
aToD.func.list = list(diag.fun, diag.fun, diag.fun, diag.fun, diag.fun, diag.fun)
aToC.fun.list =  list(aToC.fun, aToC.fun, aToC.fun, aToC.fun, aToC.fun, aToC.fun)
diag.func.list =  list(diag.fun, diag.fun, diag.fun, diag.fun, diag.fun, diag.fun)
treat.func.list = list(treat.fun,treat.fun, treat.fun, treat.fun, treat.fun,treat.fun)
svr.func.list = list(svr.fun, svr.fun, svr.fun, svr.fun, svr.fun, svr.fun)
failed_svr.func.list = list(failed_svr.fun, failed_svr.fun, failed_svr.fun, failed_svr.fun, failed_svr.fun, failed_svr.fun)
## order of transitions between blocks: 1-2(chronic), 1-3(acute_diagnosis), 2-3 (diag), 3-4 (trat), 1-6 (sr), 2-6 (sr), 3-6 (sr), 4-6 (svr), 4-5 (failed), 5-6 (svr)  
tran.blocks.func = list(aToC.fun.list, aToD.func.list, diag.func.list, treat.func.list, failed_svr.func.list, sr.func.list, sr.func.list, sr.func.list,svr.func.list, svr.func.list)

#############################################
## Setting up the transitions between blocks
# i = 5: Block 1--> 6, i = 6: 2--> 6, i = 7: 3--> 6, i = 8: 4--> 6, i = 9: 5--> 6
#############################################
## creating a set of positions
poss =  which(block.tran.mat == 1, arr.ind =  TRUE)
for (i in 1:sum(block.tran.mat)) #For each transition between blocks
{
  inter = intersect(block[[poss[i,1]]][2,], block[[poss[i,2]]][2,]) #between fibrosis levels in the two blocks taking part in i-transition
  
  pos.start   = subset(block[[poss[i,1]]][1,], block[[poss[i,1]]][2,]>= inter[1])
  pos.end     = subset(block[[poss[i,2]]][1,], block[[poss[i,2]]][2,]>= inter[1])
  
  ## Setting the transition functions to be the ones defined for in the param file
  for (j in 1:length(pos.start)) # j go through the stage vertically
  {
    tf[[pos.start[j], pos.end[j]]] = tran.blocks.func[[i]][[j]]
  }
}
########################################################
# tran.End.Stages.func = list(aToC.fun.list, aToD.func.list, diag.func.list, treat.func.list, failed_svr.func.list, sr.func.list, sr.func.list, sr.func.list,svr.func.list, svr.func.list)
# index_end0 = c(31,31, 32, 33, 34,31, 32,33, 34, 35)
# index_end1 = c(32,33, 33, 34, 35,36, 36, 36, 36, 36)
# index_end.stages = rbind(cbind(index_end0, index_end1),cbind(index_end0 + 6 * 1, index_end1 + 6 * 1),cbind(index_end0 + 6 * 2, index_end1 + 6 * 2) )
# 
# for(i in 1:sum(block.tran.mat) ) #For each transition between blocks
# tf[[index_end.stages[[i,1]] ,index_end.stages[[i,2]]]] = tran.End.Stages.func[[i]][[1]]
# 
# while(i <= 2 * sum(block.tran.mat) - 1) #For each transition between blocks
# {
#   i = i + 1
#   j = i - 10
#   tf[[index_end.stages[[i,1]] ,index_end.stages[[i,2]]]] = tran.End.Stages.func[[j]][[1]]
# }
# while(i <= 3 * sum(block.tran.mat) - 1) #For each transition between blocks
# {
#   i = i + 1
#   j = i - 20
#   tf[[index_end.stages[[i ,1]] ,index_end.stages[[i,2]]]] = tran.End.Stages.func[[j]][[1]]
# }

############################################
# for (i in 0:2)
# {
#   tf[[31 + 6 * i, 32 + 6 * i]] = aToC.fun 
#  # tf[[31 + 6 * i, 33 + 6 * i]] = diag.fun 
#   #tf[[32 + 6 * i, 33 + 6 * i]] = diag.fun
#   # tf[[33 + 6 * i, 36 + 6 * i]] = sr.fun
#   # tf[[31 + 6 * i, 36 + 6 * i]] = sr.fun
#   # tf[[32 + 6 * i, 36 + 6 * i]] = sr.fun
#   tf[[35 + 6 * i, 36 + 6 * i]] = svr.fun
#   tf[[33 + 6 * i, 34 + 6 * i]] = treat.fun
#   tf[[34 + 6 * i, 36 + 6 * i]] = svr.fun
#   tf[[34 + 6 * i, 35 + 6 * i]] = failed_svr.fun
# }
######################################################
#spontaneous recovery in DC, HCC and LT
for (i in 31:33)
  tf[[i, 36]] = sr.fun
for (i in 37:41)
  tf[[i, 42]] = sr.fun
for (i in 43:47)
  tf[[i, 48]] = sr.fun
########################################################
# Acute to chronic in DC, HCC, and LT
for (i in 0:2)
  tf[[31 + 6 * i, 32 + 6 * i]] = aToC.fun
########################################################
########################################################
# Chronic to Diagnosed
for (i in 0:2)
  tf[[32 + 6 * i, 33 + 6 * i]] = diag.fun
#######################################################s
# Diagnosed to treatment
for (i in 0:2)
  tf[[33 + 6 * i, 34 + 6 * i]] = treat.fun

for (i in 0:2)
{
  tf[[34 + 6 * i, 36 + 6 * i]] = svr.fun
  tf[[35 + 6 * i, 36 + 6 * i]] = svr.fun
}
for (i in 0:2)
  tf[[34 + 6 * i, 35 + 6 * i]] = failed_svr.fun

#######################################################
## mortality is from any state except the last one
########################################################

source("auto_mort.R")




