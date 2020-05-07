print("Please enter the svr scenario: 1. scenario_svr, 2. scenari_svr_noPl, 3.scenario_svr_fig3, 4. scenario_allcured")
scenario_svr <- readline(prompt = "Parameter_Type:")
if(scenario_svr == 1)
{
  ################################################################
  ###### diag: yearly Atb screening
  ################################################################
  svr.fun<- function(bl,history)
  {
    # identify the state we are in:
    index.fun = stage_identifier(t, lr, bl, history)
    index_block=index.fun[1]
    index_fs = index.fun[2]
    index_state = index.fun[3]
    index_previous_state = index.fun[4]
    ###########################################################
    if(bl["Genotype"] == 0) # G1
    { 
      if(index_fs!= 5) # non cirrhotics
      {
        p = (0.98) # probability of SVR
        w=runif(1,0,1)
        ifelse(w <= p && index_previous_state != (index_state - 1), 12/52, 999) 
        
      }
      else #if cirrhotics
      {
        p = (0.92)  # probability of SVR
        w = runif(1,0,1)
        ifelse( w <= p && index_previous_state != (index_state - 1), 12/52, 999)
      }
    }   
    ################################
    
    else if(bl["Genotype"] == 1) # G2
    { 
      if(index_fs != 5) # non cirrhotics
      {
        p = (0.937)
        w = runif(1,0,1)
        ifelse(w <= p && index_previous_state != (index_state - 1), 12/52, 999)
      }
      else #if cirrhotics
      {
        p=(0.818)
        w=runif(1,0,1)
        ifelse(w <= p && index_previous_state != (index_state - 1), 12/52, 999)
      }
    }
    
    ############################## 
    else if(bl["Genotype"] == 2) # G3
    { 
      if(index_fs != 5) # non cirrhotics
      {
        p = (0.91)
        w=runif(1,0,1)
        ifelse(w <= p && index_previous_state != (index_state - 1), 12/52, 999)
      }
      else #if cirrhotics
      {
        p = (0.68)
        w=runif(1,0,1)
        ifelse(w <= p && index_previous_state != (index_state -1 ), 12/52, 999)
      }
    }
    ################################# 
    else if(bl["Genotype"] == 3) # G4
    { 
      if(index_fs != 5) # non cirrhotics
      {
        p = (0.95)
        w = runif(1,0,1)
        ifelse(w <= p && index_previous_state != (index_state - 1), 12/52, 999)
      }
      else #if cirrhotics
      {
        p = (0.61)
        w=runif(1,0,1)
        ifelse(w <= p && index_previous_state != (index_state - 1), 12/52, 999)
      }
    }          
    
  }
  # 
  # 
  # svr.fun<-function(bl, history){
  #   svr_ch<-xx
  #   p=svr_ch
  #   w=runif(1,0,1)
  #   if(w<=p){
  #     if(sum(history)<=1) (12/52+0.5)
  #     else 
  #     {if (bl["Genotype"]==0 | bl["Genotype"]==1) (48/52+0.5)
  #      else (24/52 +0.5)
  #     }
  #   }
  #   
  #   else 999
  # }
}
####################################################################################
if(scenario_svr == 2)
{
  svr.fun<- function(bl,history)
  {
    # identify the state we are in:
    index.fun = stage_identifier(t, lr, bl, history)
    index_block=index.fun[1]
    index_fs = index.fun[2]
    index_state = index.fun[3]
    index_previous_state = index.fun[4]
    #######################################
    PI<-2.17
    ff<-0  
    ifelse(index_fs==1 | index_fs==2 |index_fs==3 |index_fs==4, ff<-1, ff<-0.74)
    
    { 
      if(bl["Genotype"]==1 | bl["Genotype"]==0) # G4
      { 
        if(sum(history)<1)
        {
          p=(-(0.8-(ff*0.26))/1*sum(history)+0.8)
          w=runif(1,0,1)
          if(w<=p) (12/52+0.5)
          else 999 
        }
        else #if (sum(history)>=1)
        {
          p=(ff*0.26+0*t)
          w=runif(1,0,1)
          if(w<=p) (48/52+0.5)
          else 999
        }
      }
      #############################################
      else #if(bl["Genotype"]==2 |bl["Genotype"]==3) 
      { 
        if(sum(history)<1)
        {
          p=(-(0.8-(ff*0.5))/1*sum(history)+0.8)
          w=runif(1,0,1)
          if(w<=p) (12/52+0.5)
          else 999 
        }
        else #sum(history>=1)
        {
          p=(ff*0.5+0*t)
          w=runif(1,0,1)
          if(w<=p)(12/52+0.5)
          else 999
        }
      }              
    }
    
  }
  
}
if(scenario_svr == 3)
{
  svr.fun[1]<- function(bl,history){
    # identify the state we are in:
    m=gems:::auxcounter(ns) # matrix of all transitions (without the terminal states)
    if(all(history==0)) {index_trans=1} else {index_trans=max(which(history>0))}
    #index_trans=74
    index_state=which(m==index_trans, arr.ind=TRUE)[2] # to get the colum number(arrival state) that corresponds to the last transition
    
    #return in which fibrosis state is a certain state (which level of FIbrosis prog?)
    ref0=cbind(rep(1:blocks_size, blocks_nb), 1:(ns-3))
    index_fs.fun = function(state, ref=ref0){
      ref[which(ref[,2] == state, arr.ind = TRUE),1]
    }
    index_fs=index_fs.fun(state = index_state, ref = ref0)
    #print(c("we are now in level/state", index_fs,index_state))
    PI<-1 #Telaprevir #2.17 for Boceprevir
    ff<-0  
    ifelse(index_fs==1 | index_fs==2 |index_fs==3 |index_fs==4, ff<-1, ff<-0.74)
    
    { if(bl["Genotype"]==1) # G4
    { if(sum(history)<1)
    {p=(-(0.8-(ff*0.26))/1*sum(history)+0.8)
    w=runif(1,0,1)
    if(w<=p) (12/52+0.5)
    else 999 
    }
      else #if (sum(history)>=1)
      {p=(ff*0.26+0*t)
      w=runif(1,0,1)
      if(w<=p) (48/52+0.5)
      else 999
      }
    }
      
      
      else if(bl["Genotype"]==0) # G1
      { if(sum(history)<1)
      {p=(-(0.8-(ff*PI*0.26))/1*sum(history)+0.8)
      w=runif(1,0,1)
      if(w<=p) (12/52+0.5)
      else 999 
      }
        else #sum(history>=1)
        {p=(ff*PI*0.26+0*t)
        w=runif(1,0,1)
        if(w<=p) (48/52+0.5)
        else 999
        }
      }
      
      else #if(bl["Genotype"]==2 |bl["Genotype"]==3) 
      { if(sum(history)<1)
      {p=(-(0.8-(ff*0.5))/1*sum(history)+0.8)
      w=runif(1,0,1)
      if(w<=p) (12/52+0.5)
      else 999 
      }
        else #sum(history>=1)
        {p=(ff*0.5+0*t)
        w=runif(1,0,1)
        if(w<=p)(12/52+0.5)
        else 999
        }
      }              
    }
    
    
    
    svr.fun[2]<- function(bl,history){
      # identify the state we are in:
      m=gems:::auxcounter(ns) # matrix of all transitions (without the terminal states)
      if(all(history==0)) {index_trans=1} else {index_trans=max(which(history>0))}
      #index_trans=74
      index_state=which(m==index_trans, arr.ind=TRUE)[2] # to get the colum number(arrival state) that corresponds to the last transition
      
      #return in which fibrosis state is a certain state (which level of FIbrosis prog?)
      ref0=cbind(rep(1:blocks_size, blocks_nb), 1:(ns-3))
      index_fs.fun = function(state, ref=ref0){
        ref[which(ref[,2] == state, arr.ind = TRUE),1]
      }
      index_fs=index_fs.fun(state = index_state, ref = ref0)
      #print(c("we are now in level/state", index_fs,index_state))
      PI<-1.64 #Telaprevir #2.17 for Boceprevir
      ff<-0  
      ifelse(index_fs==1 | index_fs==2 |index_fs==3 |index_fs==4, ff<-1, ff<-0.74)
      
      { if(bl["Genotype"]==1) # G4
      { if(sum(history)<1)
      {p=(-(0.8-(ff*0.26))/1*sum(history)+0.8)
      w=runif(1,0,1)
      if(w<=p) (12/52+0.5)
      else 999 
      }
        else #if (sum(history)>=1)
        {p=(ff*0.26+0*t)
        w=runif(1,0,1)
        if(w<=p) (48/52+0.5)
        else 999
        }
      }
        
        
        else if(bl["Genotype"]==0) # G1
        { if(sum(history)<1)
        {p=(-(0.8-(ff*PI*0.26))/1*sum(history)+0.8)
        w=runif(1,0,1)
        if(w<=p) (12/52+0.5)
        else 999 
        }
          else #sum(history>=1)
          {p=(ff*PI*0.26+0*t)
          w=runif(1,0,1)
          if(w<=p) (48/52+0.5)
          else 999
          }
        }
        
        else #if(bl["Genotype"]==2 |bl["Genotype"]==3) 
        { if(sum(history)<1)
        {p=(-(0.8-(ff*0.5))/1*sum(history)+0.8)
        w=runif(1,0,1)
        if(w<=p) (12/52+0.5)
        else 999 
        }
          else #sum(history>=1)
          {p=(ff*0.5+0*t)
          w=runif(1,0,1)
          if(w<=p)(12/52+0.5)
          else 999
          }
        }              
      }
      
    }
    
  }
  
  
  
  svr.fun<-function(bl, history){
    svr<-xx
    w=runif(1,0,1)
    if(w<=svr) 4/52
    else 999
  }
}

if(scenario_svr == 4)
{
  ################################################################
  ### SVR rate = 100%
  ##################################################################
  
  svr.fun<-function(bl, history){
    
    # identify the state we are in:
    m=gems:::auxcounter(ns) # matrix of all transitions (without the terminal states)
    if(all(history==0)) {index_trans=1} else {index_trans=max(which(history>0))}
    #index_trans=74
    index_state=which(m==index_trans, arr.ind=TRUE)[2] # to get the colum number(arrival state) that corresponds to the last transition
    index_previous_state=which(m==index_trans, arr.ind=TRUE)[1] # to get the row number(departure state) that corresponds to the last transition
    
    svr<-1
    w=runif(1,0,1)
    ifelse(w<=svr && index_previous_state != (index_state -1),4/52, 999 ) 
    
    
  }
  
  
}