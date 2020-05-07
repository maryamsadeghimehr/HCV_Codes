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
cat('\n \t Version: 10')
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
#setwd("~/HCV_project_2017/hcv_git/hcv_modeling")
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
start.time <- Sys.time()
current_dir = getwd()
########################################################
cohortSize <- 50 # specify how many patients to simulate
maxTime <- 100 # specify duration of simulation in years
#ans <- readline(prompt="Please enter the treatment scenario you are interested to study: Treating patients in (and after) stage F0, F1, F2, F3, F4 or BL (baseline) \t")
ans <- "Fg"
diag_scenario = c("Baseline","Birth","IDU", "ex_IDU")
diag_ind = 1
###########################################################
data <- read.table("baseline_data.txt", header = T)
bm <- read.table("Mort_CH_Swiss.txt", header=T) 

all_cases = matrix(NA,1,52 + 13)
for (case in 1:nrow(data))
{
  print(case)
  #####################Source Codes ###########################
  source("baseline.R") # call file specifying baseline characteristics
  source(sprintf("scenario_treat%s.R",ans))
  source("hazard_fun.R") # call file specifying the hazard functions
  source("params_cov.R") #choose parameters file
  source("TransitionTime.R")
  source("sim_cohort.R") # call the cohort simulator
  #########################################################
  mainDir <- wd
  today <- Sys.Date()
  dd = format(today, format = "%B%d%Y")
  subDir <- sprintf("outputDirectoryFgscenario_v1_%s_%s",dd,ans)
  dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
  wd = setwd(file.path(mainDir, subDir))
  mainDir <- wd
  setwd(file.path(mainDir, subDir))
  # print(getwd())
  write.table(cohort@time.to.state, sprintf("treatFg_v1%s_%s",case, dd))
  ########################################################
  
  write.table(cbind(cohort@time.to.state, bl), sprintf("treatFg_baseline_v1%s_%s",case, dd)) 
  all_cases = rbind(as.matrix(all_cases),as.matrix(read.table(sprintf("treatFg_baseline_v1%s_%s",case, dd))))
  setwd(current_dir)
}


#####################################

end.time <- Sys.time()
cat('***********************************************************\n')
print("Required time: ")
print(end.time - start.time)
cat('***********************************************************\n')
today <- Sys.Date()
dd = format(today, format = "%B%d%Y")

write.table(all_cases, sprintf("All_cases_Fg_v1%s_%s",dd, ans))
write.table(all_cases, sprintf("All_cases_Fg_v1%s_%s.csv",dd, ans))

