require(copula)

source("CalCopula.R")

mycsurv<- function(p1,p2,copulafam,copulaparam){
  out<- dev_Copula(copulafam = copulafam,param = copulaparam,p1 =p1,p2 =p2,mode = "2")
  return(out)
  
} 
simSCRmul<- function (n  = 100,cens.rate = 10,K=7,
                      copulafam="frank",tau.alpha=0.5,
                      tau.theta=seq(0.85,0.2,length.out=K), 
                     params=list(rate1=1, rate2=0.6)
) 
{
  
  alpha<- Calitau(copulafam,tau.alpha)
  theta<- sapply(tau.theta,Calitau,copulafam=copulafam)
  death<- rexp(n =n,rate = params$rate2)
  surv.death<- pexp(death,rate = params$rate2,lower.tail = FALSE)
  
  copula_type <- archmCopula(family= copulafam,param = alpha,dim = K)
  U <- copula::rCopula(n, copula_type) #Generate uniform random numbers from the copula
  
  # Transform uniform random numbers to samples from the desired distribution
  # Apply inverse transform method using marginal CDFs
  intermediate <- u.inter<- matrix(NA, nrow = n, ncol = K)
  for (k in 1:K) {
    for (i in 1:n) {
      u.inter[i,k]<- uniroot(f = function(xx){
        mycsurv(p1=xx,p2=surv.death[i],copulafam=copulafam,copulaparam=theta[k])-
          U[i, k]}, interval = c(0, 1))$root
    }
    
    
    
    intermediate[, k] <- qexp(u.inter[,k],rate = params$rate1,lower.tail = FALSE) #  inverse of marginal CDF
  }
  colnames(intermediate)<- paste0("T",1:k)
  # GGally::ggpairs(data.frame(cbind(intermediate,death)))
  # psych::pairs.panels(data.frame(cbind(intermediate,death)),
  #                     smooth = FALSE,scale = FALSE,density = TRUE,
  #                     ellipses = FALSE,lm=FALSE,method="kendall",
  #                     cor=TRUE,cex.cor = 1.1,stars=TRUE,
  #                     rug = FALSE,
  #                     hist.col = "lightblue1",hist.border = "black")
  
  # M<- cor(data.frame(cbind(intermediate,death)),method = "kendall")
  # pmat<- corrplot::cor.mtest(data.frame(cbind(intermediate,death)),method = "kendall")$p
  # 
  # # library(NatParksPalettes)
  # # col<-natparks.pals( "Olympic",n=6,type = "continuous",direction = -1)
  # 
  # library(RColorBrewer)
  # col<- brewer.pal(6,"BrBG")
  # corrplot::corrplot(corr = M,type="upper",
  #                    method="square",diag=FALSE,
  #                    col.lim = c(0.15,0.85),is.corr = FALSE,
  #                    p.mat =pmat,col = col)
                   
  
  C <- runif(n,min = 0 ,max = cens.rate)
  
  T2obs<- pmin(death,C)
  T1obs<- apply(intermediate,2,function(xx)pmin(xx,T2obs))
  delta<- 1*(T1obs==intermediate)
  delta.tilde<- 1*(death< C)
  
  # cat("observed rate:", apply(delta,2,sum)/n,"for intermediate events,",
  #     sum(delta.tilde)/n, "for T2","\n")

  
  data<- list(T1= T1obs,T2= T2obs,event1= delta, event2= delta.tilde,T2true= death)
  return(data)
  
}
