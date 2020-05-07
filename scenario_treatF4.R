### TREAT ONLY F4

#########treatment uptake 
treat.fun <- function(bl, history){
  
  #identify the state we are in:
  m=gems:::auxcounter(ns) # matrix of all transitions (without the terminal states)
  if(all(history==0)) {index_trans=1} else {index_trans=max(which(history>0))}
  index_state=which(m==index_trans, arr.ind=TRUE)[2] # to get the colum number(arrival state) that corresponds to the last transition
  
  #return in which fibrosis state is a certain state (which level of FIbrosis prog?)
  ref1 = cbind(rep(1:blocks_size, blocks_nb), 1:(ns - 22))
  ref2 = cbind(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6), c(31,37,43,32,38,44,33,39,45,34,40,46,35,41,47,36,42,48))
  ref0 = rbind(ref1,ref2)
  index_fs.fun = function(state, ref=ref0){
    ref[which(ref[,2] == state, arr.ind = TRUE),1]
  }
  index_fs=index_fs.fun(state = index_state, ref = ref0)
  # print(c("we are now in level/state", index_fs,index_state))
  
  if(index_fs>=5&& (bl["BirthYear"] + bl["Age"] + sum(history) + t > 2014)) 24/52 # treat in F4
  
  else 666
  
}