require(copula)

simSCRnp<- function (n  = 100,cens.rate = 5,copulafam="frank",tau, 
                     params=list(marginsDist = rep("exp",2),
                                 rate1=1, 
                                 rate2=1)
) 
{
  
  if(copulafam=="frank")alpha<- iTau(frankCopula(),tau)
  if(copulafam=="gumbel")alpha<- iTau(gumbelCopula(),tau)
  if(copulafam=="clayton")alpha<-iTau(claytonCopula(),tau)
  if(copulafam=="joe")alpha<- copula::iTau(copula::joeCopula(),tau)
  if(copulafam=="amh")alpha<- copula::iTau(copula::amhCopula(),tau)
  
  biSurvTimes<-  rMvdc(n,mvdc(copula =archmCopula(family= copulafam,
                                                  param = alpha,
                                                  dim = 2), 
                              margins = params$marginsDist, 
                              paramMargins = c(list(list(rate = params$rate1)),
                                               list(list(rate = params$rate2)))))
  T1<- biSurvTimes[,1]
  T2<- biSurvTimes[,2]
  
  C <- runif(n,min = 0 ,max = cens.rate)
  
  T2obs<- pmin(T2,C)
  T1obs<- pmin(T1,T2obs)
  delta1<- 1*(T1<= T2obs)
  delta2<- 1*(T2<= C)
  
  cat("observed rate:", sum(delta1)/n,"for T1,",sum(delta2)/n, "for T2","\n")
  
  
  data<- data.frame(T1= T1obs,T2= T2obs,event1= delta1, event2= delta2)
  return(data)
  
}
