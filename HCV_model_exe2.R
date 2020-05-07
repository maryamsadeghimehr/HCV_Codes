#############Executable file ######################################
#########################################
### Scenario: We need to choose the treatment scenario that we want to apply 
### Simulation starts at time of HCV infection
### This file include all the required packages, and call all the required sources
### In this file the Cohort Size and the duration of simulation are defined. 
#####################################################################
cat("\014")  
cat('\n')
cat('***********************************************************\n')
cat('Mathematical Modeling of Hepatitis C Virus Progression Disease\n')
cat('\t \t 2017')
cat('\n')
cat ('\n \t Version: 10')
cat('\n Requires files: \n')
cat('simulate_cohort.r: To do the simulation \n')
cat('baseline.r: Identify the baseline chracteristics \n')
cat('hazard_fun.r: Calls file specifying the hazard functions \n')
cat('params_cov.r: To choose parameters file \n')
cat('scenario_treatFx.r: To choose the treatment scenario \n')
cat('simulate_cohort.r: To do the simulation \n')
cat('TransitionTime.r: Shows all the possible transitions\n')
cat('***********************************************************\n')
#############################################################
setwd("~/HCV_project_2017/hcv_git/hcv_modeling")
rm(list =  ls())
wd =  getwd()
######################### Required Packages #################
require(MASS)
require(msm)
require(mstate)
require(survival)
require(splines)
require(plyr)
require(gems)
#########################################################
start.time<-Sys.time()
########################################################
cohortSize <- 50 # specify how many patients to simulate
maxTime<-100 # specify duration of simulation in years
setwd("../Required_Data")
data <- read.table("baseline_data.txt", header=T)
setwd("~/HCV_project_2017/hcv_git/hcv_modeling")
all_cases = matrix(NA,1,52)
setwd("~/HCV_project_2017/hcv_git/hcv_modeling")
#for (case in 1: nrow(data))
  for (case in 219: nrow(data))
{
  print(case)
    #####################Source Codes ###########################
    source("baseline.r") # call file specifying baseline characteristics
    setwd("~/HCV_project_2017/hcv_git/hcv_modeling")
    
    cat('***********************************************************\n')
    #ans <- readline(prompt="Please enter the treatment scenario you are interested to study: Treating patients in (and after) stage F0, F1, F2, F3, F4 or BL (baseline) \t")
    ans <- "F2"
    source(sprintf("scenario_treat%s.r",ans))
    #source("scenario_treatF3.r")
    source("hazard_fun.r") # call file specifying the hazard functions
    source("params_cov.r") #choose parameters file
    source("TransitionTime.r")
    source("sim_cohort.r") # call the cohort simulator
  cat('***********************************************************\n')
  ######################################################
  setwd(wd)
  setwd("../Tables")
  #########################################################
  ########################################################
  write.table(cohort@time.to.state, sprintf("SS_%s",case)) 
  all_cases = rbind(as.matrix(all_cases),as.matrix(read.table(sprintf("SS_%s",case))))
  setwd("~/HCV_project_2017/hcv_git/hcv_modeling")
  
}
#####################################

end.time <- Sys.time()
cat('***********************************************************\n')
print("Required time: ")
print (end.time - start.time)
cat('***********************************************************\n')
cat('You can find all the required data in the Table "all_cases"!! \n')
answer <- readline(prompt = "Do you want to save all_cases in Required_Data folder?")
if (answer == "y")
{
  setwd("../Required_Data")
  write.table(all_cases, "All_cases")
}
setwd("~/HCV_project_2017/hcv_git/hcv_modeling")
cat('***********************************************************\n')
tt <- readline(prompt="Do you want to save the results as a table for post processing (y/n)?")
if (tt == "y")
{
  setwd("../Outcomes")
  ans <- readline(prompt="Please Enter a File Name:")
  write.table(cohort@time.to.state,ans) # save the cohort
  Total = cbind(cohort@time.to.state, bl)
  write.table(Total, sprintf("Total_results_%s.txt",ans))
}
print("Do you want to do the post processing?")
print("Do you want to do the post processing?")
answer <- readline(prompt="Do you want to do the post processing(y/n)?")
if (answer == "y")
{

setwd("../Outcomes")
source("outcomes.r")
source("out_bc.R")
source("barPlot.r")
source("out_baseline.r")
  
}

