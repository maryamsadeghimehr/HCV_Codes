Treat_scenario <- 1
if (Treat_scenario == 1)
{
#########treatment uptake 
# In this scenario individuals will be treated if they are in 2014 or more, and they have fibrosis F2 or more, 
# also if they have fibrosis less than F2, they will be treated if they are in 2018 or more.
treat.fun <- function(bl, history)
{
  source("index_block_fs_stage.r")
  lr = 0
  index = stage_identifier(lr , bl, history)
  index_fs = index[2]
  state = index[3]
  
  if (is.element(state, c(31:36)) )
    index_fs = 6 # DC
  
  if (is.element(state, c(37:42)) )
    index_fs = 7 # HCC
  
  if (is.element(state, c(43:48)) )
    index_fs = 8 #LT
  #############################################
  Year_infection = bl[ "BirthYear"] + bl["Age"]
  epsi = runif(1, 0, 1)
  res <- ifelse(index_fs >= 3, (2014 - (sum(history) + Year_infection)), 
                (2018 - (sum(history) + Year_infection)))  # treat in F2 or above 
  epsi = ifelse(sum(history) + Year_infection < 2014, runif(1,0,5) , epsi)

  res <- max(0, res) + epsi
  return(res)
}
}



if (Treat_scenario == 2)
{

#########treatment uptake 
# In this scenario individuals will be treated if they are in 2014 or more, and they have fibrosis F2 or more, 
# also if they have fibrosis less than F2, they will be treated if they are in 2018 or more.
treat.fun <- function(bl, history)
{
  source("index_block_fs_stage.r")
  lr = 0
  index = stage_identifier(lr , bl, history)
  index_fs = index[2]
  state = index[3]
  
  if (is.element(state, c(31:36)) )
    index_fs = 6 # DC
  
  if (is.element(state, c(37:42)) )
    index_fs = 7 # HCC
  
  if (is.element(state, c(43:48)) )
    index_fs = 8 #LT
  #############################################
  Year_infection = bl[ "BirthYear"] + bl["Age"]
 # epsi = runif(1, 0, 1)
  epsi = ((sum(history) + Year_infection) - 2014) / (1 + (sum(history) + Year_infection - 2014))
  res <- ifelse(index_fs >= 3, (2014 - (sum(history) + Year_infection)), 
                (2018 - (sum(history) + Year_infection)))  # treat in F2 or above 
  epsi = ifelse(sum(history) + Year_infection < 2014, runif(1,0,5) , epsi)
  
  res <- max(0, res) + epsi
  return(res)
}

}

