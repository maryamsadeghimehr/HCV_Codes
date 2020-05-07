#################################################
###
###
################################################
#treatment uptake 
treat.fun <- function(bl, history)
  {
  source("index_block_fs_stage.r") # it identify that in which block and stage we are in
  index.fun = stage_identifier(lr, bl, history)
  index_fs = index.fun[2]
  if(index_fs>= 4 && (bl["BirthYear"] + bl["Age"] + sum(history) + t > 2014)) res <- rep(24/52, length(t)) # treat in F3
  
  else res <- rep(666, length(t))
  
}
