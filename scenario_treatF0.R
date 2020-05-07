#########treatment uptake 
treat.fun <- function(bl, history, t)
  {
if(bl["BirthYear"] + bl["Age"] + sum(history) + t > 2014)
  res <- rep((6 * 4)/52, length(t)) # treat in F0
  else
    res <- rep(666, length(t))
 }