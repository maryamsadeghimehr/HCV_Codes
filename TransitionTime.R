#####################################################################
### specifies which transitions are "time to transition" 
#####################################################################
#template for time to transition logic
timeToTrans <- matrix(data = FALSE, nrow = statesNum, ncol = statesNum)
####################################################################
# specifies which transitions are "time to transition" 
for (i in 1:5)
{
  timeToTrans [i, i + 5] = TRUE # progress horizantally
}

for(i in 1:49)
{
  
  if (i > 15 && i < 26)
  {
    if (i %% blocks_size == 0)
      timeToTrans [i, blocks_size * blocks_nb] = TRUE # spontaneous recovery for state F4 at any block
    else
      timeToTrans [i, blocks_size * (blocks_size - floor(i/ blocks_size)) + i] = TRUE # Spontaneous recovery
    
  }# end of if (i > 15 && i < 26)
  ########################################
  if (i < 16) 
  {
    if (i %% blocks_size == 0)
      timeToTrans [i, blocks_size * blocks_nb] = TRUE # spontaneous recovery for state F4 at any block
    else
    {
      timeToTrans [i, blocks_size * (blocks_size - floor(i/ blocks_size)) + i] = TRUE # Spontaneous recovery
    }
  }
  
}

for (i in 11:25)
  timeToTrans [i, i + 5] = TRUE
##########################################
#spontaneous recovery in DC, HCC and LT
for (i in 31:35)
  timeToTrans [i, 36] = TRUE
for (i in 37:41)
  timeToTrans [i, 42] = TRUE
for (i in 43: 47)
  timeToTrans [i, 48] = TRUE

############################################
for (i in 0:2)
{
  timeToTrans [35 + 6 * i, 36 + 6 * i] = TRUE #SVR
  timeToTrans [33 + 6 * i, 34 + 6 * i] = TRUE # Scenario.treatment
  timeToTrans [34 + 6 * i, 36 + 6 * i] = TRUE # SVR
  timeToTrans [34 + 6 * i, 35 + 6 * i] = TRUE # second treatment
}
#which(timeToTrans == TRUE, arr.ind = TRUE)

####################################################
timeToTrans[31, 32] = TRUE 
timeToTrans[37,38] = TRUE 
timeToTrans[43,44] = TRUE 
