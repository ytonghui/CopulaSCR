scrassonp<- function(data, t1.formula, t2.formula,weights, na.action,
                   equalweight = FALSE,copulafam=c("clayton","frank","joe","gumbel","amh"),
                   se=FALSE,se.method = c("bootstrap","resampling"),B=100,model=FALSE,
                   a=0,b=0,tol=1e-5,seed=NULL,positive=TRUE){
  
  Call <- match.call()
  
  se.method<- match.arg(se.method)
  copulafam<- match.arg(copulafam)
  
  indx <- match(c('data', 't1.formula', 't2.formula', 'weights,' ,'na.action'), names(Call), nomatch=0)
  if (indx[1]==0) stop("a data argument is required")
  if (indx[2]==0) stop("a t1.formula argument is required")
  if (indx[3]==0) stop("a t2.formula argument is required")
  
  
  if (!(is.matrix(data) | is.data.frame(data))) {
    stop("data object must be a matrix or a data.frame")
  }
  
  if (indx[4]==0&indx[5]==0){
    t1data <- prepData(formula=t1.formula, data)
    t2data <- prepData(formula=t2.formula, data)
  }
  if(indx[4]==0&indx[5]!=0){
    t1data <- prepData(formula=t1.formula, data,na.action = na.action)
    t2data <- prepData(formula=t2.formula, data,na.action = na.action)
  }
  
  if(indx[4]!=0&indx[5]==0){
    t1data <- prepData(formula=t1.formula, data,weights = weights)
    t2data <- prepData(formula=t2.formula, data,weights = weights)
  }
  if(indx[4]!=0&indx[5]!=0){
    t1data <- prepData(formula=t1.formula, data,weights = weights,na.action = na.action)
    t2data <- prepData(formula=t2.formula, data,weights = weights,na.action = na.action)
  }
  
  casewt<- t1data$weights
  time1<- t1data$time
  event1<- t1data$event
  time2<- t2data$time
  event2<- t2data$event
  offset1<- t1data$offset
  offset2<- t2data$offset
  
  n<- length(time1)
  
  if (is.null(casewt)) {
    casewt <- rep(1.0, n)
  }  else {
    if (!is.numeric(casewt)) stop("weights must be numeric")
    if (any(!is.finite(casewt))) stop("weights must be finite") 
    if (any(casewt <0)) stop("weights must be non-negative")
    casewt <- as.numeric(casewt)  # transform integer to numeric
  }
  
  
  copulaparam<- assonpfit( time1 = time1,event1 = event1,time2 = time2,event2 = event2,
                           copulafam = copulafam,equalweight=equalweight, 
                           a=a,b= b,tol=tol,zeta=NULL,positive=positive)
  
  if(copulafam=="frank")taufn<- function(copulaparam) copula::tau(copula::frankCopula(copulaparam))
  if(copulafam=="gumbel")taufn<- function(copulaparam)copula::tau(copula::gumbelCopula(copulaparam))
  if(copulafam=="clayton")taufn<- function(copulaparam)copula::tau(copula::claytonCopula(copulaparam))
  if(copulafam=="joe")taufn<- function(copulaparam)copula::tau(copula::joeCopula(copulaparam))
  if(copulafam=="amh")taufn<- function(copulaparam)copula::tau(copula::amhCopula(copulaparam))
  
  tau.est<- taufn(copulaparam)
  
  out<- list(Call=Call,tau=tau.est,copulaparam=copulaparam,copulafam = copulafam)
  out$t1.formula<- t1.formula
  out$t2.formula<- t2.formula
  
  if(isTRUE(model)){
    out$data<- data
  }
  
  out$se<- se
  if(isTRUE(se)&!is.null(seed)&is.numeric(seed)){set.seed(seed)}
  
  if(isTRUE(se)){
    
    out$se.method<- se.method
    copulaparam.boot<- c()
    cat("Please wait for a while...")
    if(se.method =="resampling"){
      zetas<- c()
      for(k in 1:B){
        zeta<-rexp(n,rate=1)
        copulaparam.boot<- c(copulaparam.boot, assonpfit( time1 = time1,event1 = event1,time2 = time2,event2 = event2,
                          copulafam = copulafam,equalweight=equalweight, 
                          a=a,b= b,tol=tol,zeta=zeta,positive=positive))
        zetas<- rbind(zetas,zeta)
        
      }
      out$zetas<- zetas
      
    }else{
      bootind<- c()
      for(k in 1:B){
        indd<- sort(sample(1:n,replace = TRUE))
        
        copulaparam.boot<- c(copulaparam.boot, assonpfit( time1 = time1[indd],event1 = event1[indd],
                                                         time2 = time2[indd],event2 = event2[indd],
                                                         copulafam = copulafam,equalweight=equalweight, 
                                                         a=a,b= b,tol=tol,zeta=NULL,positive=positive))
        
        bootind<- rbind(bootind,indd)
      }
      out$bootind<- bootind
    }
    out$param.boot<- copulaparam.boot
    out$param.se<- sd(copulaparam.boot)
    out$B<- B
    out$tau.boot<- sapply(copulaparam.boot, taufn)
    out$tau.se<- sd(out$tau.boot)
  }
  class(out)<- "scrassonp"
  return(out)
  
}


assonpfit<- function( time1,event1,time2,event2,equalweight,copulafam, a,b,tol=1e-5,zeta=NULL,positive=TRUE){
  n<- length(event2)
  if(is.null(zeta)){
    zeta2<- 1
    }else{
      zeta2<- combn( zeta,m=2)
      zeta2<- zeta2[1,]*zeta2[2,]
    }
  
  t1.pair<- combn(time1,m=2)
  t2.pair<- combn(time2,m=2)
  event1.pair<- combn(event1,m=2)
  event2.pair<- combn(event2,m=2)
  
  t1.pairmin<- apply(t1.pair,2,min)
  t2.pairmin<- apply(t2.pair,2,min)
  
  pair<- ((t1.pair[1,]>t1.pair[2,])*event1.pair[2,]+
            (t1.pair[1,]<t1.pair[2,])*event1.pair[1,])
  pair<- pair*((t2.pair[1,]<t2.pair[2,]) *event2.pair[1,]
               +(t2.pair[1,]>t2.pair[2,])*event2.pair[2,])
  
  
  Delta<- (t1.pair[1,]-t1.pair[2,])*(t2.pair[1,]-t2.pair[2,])
  Delta<- 1* (Delta>0 )
  
  
  if(equalweight){
    w<- 1
  }else{
    if(is.null(a)){a<- quantile(time1,0.95)}
    if(is.null(b)){b<- quantile(time2,0.95)}
    w_inv<- sapply(1:length(t1.pairmin),function(kk)
      mean((time1>=min(a,t1.pairmin[kk]))*(time2>=min(b,t2.pairmin[kk]))))
    w<- 1/w_inv
    rm(w_inv)
  }
  
  w<- w*zeta2
  
  if(copulafam=="clayton"){
    numer<- sum(w*Delta *pair,na.rm = TRUE)
    denom<- sum(w*(1-Delta) *pair,na.rm = TRUE)
    copulaparam<- numer/denom -1
    
  }else{
    js0<- sapply(1:length(t1.pairmin),function(kk)
      mean((time1>t1.pairmin[kk])*(time2>t2.pairmin[kk])))
    if(sum(1-event2)/n>0.05){
      Gfit<- survival::survfit(survival::Surv(time2,1-event2)~1,weights = zeta)
      Gest<- approx(c(0,Gfit$time),c(1,Gfit$surv),yright = min(Gfit$surv),ties = mean,
                    xout =  t2.pairmin,method="constant",f = 0)$y
      # Gfit<- prodlim::prodlim(Surv(time2,event2)~1,reverse = TRUE,caseweights = zeta)
      # Gest<-  predict(Gfit,times = t2.pairmin)
      rm(Gfit)
      
      js<- js0/ Gest
      js[js0==0& Gest==0]<- 0
    }else{
      js<- js0
    }
    
    if(isTRUE(positive)){
      paramin<-copulabd_pos(copulafam)
    }else{
      paramin<- copulabd(copulafam)
    }
    
    
    # copulaparam<- stats::uniroot(f = function(par){
    #   scrnpGEE(gamma=par,copulafam = copulafam,w = w,js = js,pair = pair,
    #            Delta = Delta)},
    #   interval = paramin,trace = FALSE,maxiter =maxiter,tol = tol)$root
    
    copulaparam<- stats::optimize(f = function(par){
      abs(scrnpGEE(gamma=par,copulafam = copulafam,w = w,js = js,pair = pair,
                   Delta = Delta))},
      interval = paramin,tol = tol)$minimum
    
    
  }
  return(copulaparam)
  
}
scrnpGEE<- function(gamma,copulafam,w,js,pair,Delta){
  
  chr<- npCHRfit(s=js,param = gamma,copulafam = copulafam)
  mean(w *pair* (Delta-  chr/(chr+1)),na.rm = TRUE)
  
  # mean((w *(Delta-  chr/(chr+1)))[pair==1],na.rm = TRUE)
}


npCHRfit<- function(s,param,copulafam){
  if(copulafam=="amh"){
    phi<- (2*param*s +1-param)/(param*s +1-param)
  }
  if(copulafam=="clayton"){
    phi<- param +1
  }
  if(copulafam=="frank"){
    phi<- (param*s)/(1-exp(-param*s))
    phi[s==0]<- 1
  }
  if(copulafam=="gumbel"){
    phi<- 1 + (1-param)/log(s)
  }
  if(copulafam=="joe"){
    phi<-  s*(param-1+(1-s)^param)/(1-s) /(1-(1-s)^param)
    phi[s==0]<- 1
  }
  return(phi)
}

copulabd_pos<- function(copulafam){
  lb<- switch(copulafam, 
              amh = -1, 
              clayton = 0, 
              frank = 0.0001, 
              gumbel = 1, 
              joe = 1)
  ub<- switch(copulafam, 
              amh = 1, 
              clayton = 1000, 
              frank = 1000, 
              gumbel = 1000, 
              joe = 1000)
  return(c(lb,ub))
}


