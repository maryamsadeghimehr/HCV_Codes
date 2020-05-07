######################################################
## Cohort simulator
# We have 52 states. For each state all trnsition specific hazard functions and their parameters need to be specified.
# Output: Time to state (entry time for each patient into each of the states.)
######################################################

cohort <- simulateCohort(
  transitionFunction = tf, #Contain hazard functions and time to events
  parameters = pm, # Contain the parameters
  cohortSize = cohortSize, # Number of patients to be simulated
  parameterCovariance = FALSE, # Covarience matrix for parameters
  sampler.steps = 300,
  report.every = 1000,
  baseline = bl, # A matrix or data frame of baseline characteristics.
  initialState = start,
  to = maxTime ,
  absorbing = 49:52,
  timeToTransition =  timeToTrans # A logical matrix, it is true for all transitions whose transition function is specified as time until transition instead of hazard function.
  
)
##################################################################
#cat('***********************************************************\n')
#print("Time to entering each state:")
cat('***********************************************************\n')
Cohort <- cohort@time.to.state
colnames(Cohort) <- c(rep("acute",5), rep("chronic_undiag",5),rep("diag",5),rep("first_treatment",5), rep("second_treatment",5),rep("cleared",5),"DC_acute", "DC_chronic_undiag", "DC_diag","DC_first_treatment", "DC_second_treatment","DC_clear",
                      "HCC_acute", "HCC_chronic_undiag", "HCC_diag","HCC_first_treatment", "HCC_second_treatment","HCC_clear",
                      "LT_acute", "LT_chronic_undiag", "LT_diag","LT_first_treatment", "LT_second_treatment","DC_clear","B_death","IDU_death", "HIV_death","Liver_death")
#print(Cohort)
#cat('***********************************************************\n')
#cat('***********************************************************\n')
#print("Calendar time")

cat('***********************************************************\n')
Cohort_calendar_time <- bl[,"Age"] + bl[,"BirthYear"] + Cohort
#print(Cohort_calendar_time)