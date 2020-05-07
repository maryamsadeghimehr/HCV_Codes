# Find an appropriate block for states: blocks_size = 5
# We have 6 DC + 6 HCC + 6 LT + 4 Death = 22 cases
stage_identifier = function(lr = 0, bl, history)
{
  statesNum = 52
  # identify the state we are in:
  m = gems:::auxcounter(statesNum) # matrix of all transitions (without the terminal states)
  if (all(history == 0)) {index_trans = 1} else {index_trans = max(which(history > 0))}
  index_state = which(m == index_trans, arr.ind = TRUE)[2] # to get the colum number(arrival state) that corresponds to the last transition
  index_previous_state = which(m == index_trans, arr.ind = TRUE)[1] # to get the row number(departure state) that corresponds to the last transition
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
  # #######################################
  # ref = ref2
  # for (state in 31:36)
  # {
  #   ref[which(ref2[,2] == state, arr.ind = TRUE),1] * 0 + 6 # DC
  # }
  #     
  # ######################################
  # 
  # for (state in 37:42)
  #   ref[which(ref2[,2] == state, arr.ind = TRUE),1] * 0 + 7 # HCC
  #   #index_fs = 7 # HCC
  # ####################################
  # for (state in 43:48)
  #   ref[which(ref2[,2] == state, arr.ind = TRUE),1] * 0 + 8 # LT
  #  # index_fs = 8 #LT
  # ####################################
  index = c(index_block, index_fs,index_state, index_previous_state)
}

