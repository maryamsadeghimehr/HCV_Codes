# Alcohol: 
# no alcohol abuse: 0, alcohol abuse: 1, 2 more than 40
# Genotype: 0,1,2,3
# HIV: 0, 1
# MSM: 0, 1
# IDU: 0, 1
# Origin:
# Swiss, Italy, Germany, Uugoslavia, Portugal, France, Spain, Austria, Georgia, Russia, others
#     1,3,4,5,6,7,8,9,10
# age:
#  [< 21), [21 - 31), [31 - 41), [41 - 51), [51 - 61), [61 - 71), [71 - 81),  [81 - ???)
# Year of birth: 

#[1937 - 1947], [1947 - 1957], [1957 - 1967], [1967 - 1977], [1977 - 1987], [1987 - 1997], [1997 - 2007], [2007-2016]
# for the origin we can have Swiss: 0, not swiss from high prevalence: 1, not swiss from low prevalence: 2, not swiss from medium prevalence
###############################################################################
setwd("~/HCV_project_2017/hcv_git//Required_Data")
###############################################################################
library(data.table)  
baseline = CJ(c(0,1), c(0,1,2), c(0,1), c(0,1), c(0,1),c(1),c(0,1,2,3),c(1,2,3,4,5, 6, 7, 8), c(1,2,3,4,5,6,7,8))
colnames(baseline)<-c("Alcohol","Gender","HIV","MSM","IDU", "Genotype" , "Origin", "Age", "BirthYear")
###############################################################################
dd <- baseline[!(baseline$"Gender" == 0 & baseline$"MSM" == 1),] # female
 ##########################################################################
 dd <- baseline[!(baseline$"Age" == 5 & baseline$"BirthYear" == 1),]
 #########################################################################
 dd <- baseline[!(baseline$"Age" == 6 & baseline$"BirthYear" == 1),]
 dd <- baseline[!(baseline$"Age" == 6 & baseline$"BirthYear" == 2),]
 dd <- baseline[!(baseline$"Age" == 6 & baseline$"BirthYear" == 3),]
 dd <- baseline[!(baseline$"Age" == 6 & baseline$"BirthYear" == 4),]
 ##########################################################################
 dd <- baseline[!(baseline$"Age" == 7 & baseline$"BirthYear" == 1),]
 dd <- baseline[!(baseline$"Age" == 7 & baseline$"BirthYear" == 2),]
 dd <- baseline[!(baseline$"Age" == 7 & baseline$"BirthYear" == 3),]
 dd <- baseline[!(baseline$"Age" == 7 & baseline$"BirthYear" == 4),]
 ##########################################################################
 dd <- baseline[!(baseline$"Age" == 8 & baseline$"BirthYear" == 1),]
 dd <- baseline[!(baseline$"Age" == 8 & baseline$"BirthYear" == 2),]
 dd <- baseline[!(baseline$"Age" == 8 & baseline$"BirthYear" == 3),]
 dd <- baseline[!(baseline$"Age" == 8 & baseline$"BirthYear" == 4),]
 dd <- baseline[!(baseline$"Age" == 8 & baseline$"BirthYear" == 5),]
 #########################################################################
 write.table(dd, "baseline_data.txt")
 setwd("~/HCV_project_2017/hcv_git/hcv_modeling")
 
 