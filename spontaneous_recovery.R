# print("Please enter the interested spontaneous recovery scenario: 1. as baseline chracteristics, 2. A logistic decrease 3. hazard function:")
# 
# sr_scenario <- readline(prompt = "sr_scenario:")
sr_scenario = 1
if (sr_scenario == 1)
{
  sr.fun = function (bl, history)
  {
    bl["sr.time"]
  }
}
##################################################
# spontaneous recovery with logistic decrease
#################################################
if (sr_scenario == 2)
{
  sr.fun = function(bl, history)
  {
    # r_sr<-0.0000000123
    # res<-r_sr+0*t
    logisticfunc <- function(min, max, infl, slo) {min +((max - min)/(1 + (sum(history) / infl) ^ slo))}
    
    p = logisticfunc(min = 0, max = 1, infl = 0.25, slo = 2.23)
    w = runif(1,0,1)
    #if(w<=p && sum(history)<=1) (0.001)
    if(w<= p) ( 0.001)
    else 999 
  } 
}
#################################################
