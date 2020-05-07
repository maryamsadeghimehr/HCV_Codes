To run the simulation:
Run "HCV_model_exe.R" in R

To see the description of the project:
See "HCV_projec.html"
 
 
 
############################################################ 
More details:
The states:
 Acute: 1:5, 31(DC), 37 (HCC), 43 (LT)
 undiag: 6:10, 32(DC), 38 (HCC), 44 (LT)
 diag: 11:15, 33(DC), 39 (HCC), 45 (LT)
 first treated: 16:20, 34(DC), 40 (HCC), 46 (LT)
 second treated: 21:25, 35(DC), 41 (HCC), 47 (LT)
 undetected: 26: 30, 36(DC), 42 (HCC), 48 (LT)
Death 49: 52
49: Background mortality , 50: IDU Death, 51: HIV Death, 52: Liver Death
 
 Baseline characteristics

 A matrix, "bl", which consists of the following characteristics of the patients as the baseline characteristics:
   "Alcohol", "Genotype", "Age", "Origin", "HIV", "MSM", "IDU", "Diag.time", "Treat.time", "sr.time".

 Liver progression 
 
 File name: hazard_fun.r

 All the hazard functions are defined in this file.
	Progression from F0 to F4: fp.fun (9 different scenario)
 
 Default:  liver disease progression through F0 to F4 depends on age at HCV infection, and viral load.
 Age at HCV infection:  rate of fibrosis progression per year based on age at HCV infection
 Viral load:  rate ratio: 0.1 undetectable HCV (under treatment or SVR) / detectable HCV VL
 
 Progression from F4 to DC: dc.fun
 Depends on the viral load
 Viral load:  rate ratio: 0.1 undetectable HCV (under treatment or SVR) / detectable HCV VL

 
 Progression from F4 to HCC: hcc.fun
 Depends on the viral load
 Viral load:  rate ratio: 0.38 undetectable HCV (under treatment or SVR) / detectable HCV VL
 rate ratio ? exp (lr) (where exp (lr) is set to 1)
 
 Progression from DC to HCC: dchcc.fun
 Fixed value exp(lr)
 
 .	Progression from DC to LT: dcLT.fun
 .	Progression from HCC to LT: hccLT.fun
 


 
 Cascade of care
 
 .	aToc.fun: Progression from acute to chronic
 fixed time to event: 6/12
 .	sr.fun: Spontaneous recovery
 -	fixed time to event (time to spontaneous clearance is calculated before the simulation and assigned as a baseline)
 -	fixed time to event (spontaneous recovery with logistic decrease)
 
 .	diag.fun: diagnosis time
 -	diagnosis time can be calculated before the simulation and assigned as a baseline
 -	fixed time to event which depends on probability of detection with the specific test we use. It also can depend on the period of screening  and the probability of screening
 .	treat.fun: How long after the HCV diagnosis treatment will be taken 
 -	treat one month after reaching F0, F1, F2, F3 or F4
 -	define as a baseline characteristic and treat everyone a particular time after the detection
 .	SVR.fun
 Time to event, probability of SVR depends on Genotype    

 Death:
 
 .	mort.fun
 -	mortality from DC/HCC after the first year
 -	mortality from DC/HCC during the first year
 -	adjustment for HIV related mortality
 -	adjustment for IDU related mortality

 mortality_adjustment.r

When we set the simulation time to 100, we are sure that all the people who die due to HCV without entering to cleared states, will die. 
 Looking at the mortality results it can be observed that the ones who do not die during the simulation are basically the ones who recover the HCV at the age less than 100,
 or the ones who are failed in the first line treatment and do not go for the second course of the treatment. 
 This file aim at the mortality of this group and try to fit the time of death to these individulas.

        