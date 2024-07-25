library(survival)
scrsurv <- function(data, fit, t1.formula, t2.formula,weights, na.action,
                    tau,copulafam=c("clayton","frank","joe","gumbel","amh"), 
                    method= c("JFKC","LRA","FJC","Ker"),surv2km=TRUE,
                    se.method = c("bootstrap","resampling"),B=100,
                    confint=TRUE, conftype=1,model=TRUE,seed=NULL ) {
  Call <- match.call()
  
  method<- match.arg(method)
  method<- substr(toupper(method),1,1)
  copulafam<- match.arg(copulafam)
  se.method<- match.arg(se.method)
  indx <- match(c('data','fit' ,'t1.formula', 't2.formula', 'weights,' ,'na.action'), names(Call), nomatch=0)
  
  if (indx[2]==0) {
    if (indx[1]==0) stop("a data argument is required")
    if (indx[3]==0) stop("a t1.formula argument is required")
    if (indx[4]==0) stop("a t2.formula argument is required")
  }
  
  
  if (indx[2]!=0) {
    if(class(fit)!= "scrassonp"){
      stop("fit shall be 'scrassonp' class.")
    }else if (is.null(fit$data)&indx[1]==0) {
      stop("a data argument is required")
    }else if (!is.null(fit$data)&indx[1]==0){
      data<- fit$data
    }
    
    if (indx[3]==0)  t1.formula<- fit$t1.formula
    if (indx[4]==0) t2.formula<- fit$t2.formula
    
    copulafam<- fit$copulafam
    tau<- fit$tau
    copulaparam<- fit$copulaparam
    
  }
  if (!(is.matrix(data) | is.data.frame(data))) {
    stop("data object must be a matrix or a data.frame")
  }
  
  if (indx[5]==0&indx[6]==0){
    t1data <- prepData(formula=t1.formula, data = data)
    t2data <- prepData(formula=t2.formula, data = data)
  }
  if (indx[5]==0&indx[6]!=0){
    t1data <- prepData(formula=t1.formula, data = data,na.action = na.action)
    t2data <- prepData(formula=t2.formula, data = data,na.action = na.action)
  }
  
  if (indx[5]!=0&indx[6]==0){
    t1data <- prepData(formula=t1.formula, data = data,weights = weights)
    t2data <- prepData(formula=t2.formula, data = data,weights = weights)
  }
  if(indx[5]!=0&indx[6]!=0){
    t1data <- prepData(formula=t1.formula, data = data,weights = weights,na.action = na.action)
    t2data <- prepData(formula=t2.formula, data = data,weights = weights,na.action = na.action)
  }
  
  if(isTRUE(confint)&!is.null(seed)&is.numeric(seed)){set.seed(seed)}
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
  
  if(method!="K"){
    if(copulafam=="frank")copulaparam<- iTau(frankCopula(),tau)
    if(copulafam=="gumbel")copulaparam<- iTau(gumbelCopula(),tau)
    if(copulafam=="clayton")copulaparam<-iTau(claytonCopula(),tau)
    if(copulafam=="joe")copulaparam<-iTau(joeCopula(),tau)
    if(copulafam=="amh")copulaparam<-iTau(amhCopula(),tau)
  }
  if(method!="J"){se.method <- "bootstrap"}
  if(se.method !="resampling"){se.method <- "bootstrap"}
  
  if(!isTRUE(confint)|conftype%in% 1:2){
    sout<- scrsurvfit(time1 = time1,event1 = event1,time2 = time2,event2 = event2,
                       copulaparam = copulaparam,copulafam = copulafam,
                       model = model,casewt=casewt,surv2km=surv2km, 
                       method=method,confint=confint, se.method=se.method,
                       B=B,conftype=conftype)
  }
  if(isTRUE(confint)&conftype%in% 3:4){
    if (!isTRUE(fit$se)){
      stop("a se argument in fit shall be TRUE")
    }
    if (fit$se.method!= se.method ){
      stop("se.method arguments in fit and scrsurv function shall be consistent.")
    }
    
    if (fit$B!= B ){
      stop("'B' in 'CopulaSCRnp' and scrsurv function shall be consistent.")
    }
    
    sout<- scrsurvfit(time1 = time1,event1 = event1,time2 = time2,event2 = event2,
                      copulaparam = copulaparam,copulafam = copulafam,
                      copulaparams = fit$param.boot,bootind = fit$bootind,zetas = fit$zetas,
                      model = model,casewt=casewt,surv2km=surv2km, 
                       method=method,confint=confint, se.method=se.method,
                       B=B,conftype=conftype)
  }
  if(method!="K"){
    sout$copulaparam<- copulaparam
    sout$copulafam<- copulafam
    sout$tau<- tau
  }
  sout$Call<- Call
  sout$t1.formula<- t1.formula
  sout$t2.formula<- t2.formula
  if(isTRUE(confint)){sout$se.method<- se.method  }
  if(isTRUE(confint)&conftype%in% 3:4){
    sout$se.method<- se.method
    sout$copulaparams = fit$param.boot
    sout$conftype<- conftype
    # if(se.method == "bootstrap") {
    #   sout$bootind = fit$bootind
    # }
    # if(se.method == "resampling") {
    #   sout$zetas = fit$zetas
    # }
    
  }
  
  if(isTRUE(model)){
    sout$t1data<- t1data
    sout$t2data<- t2data
  }
  class(sout)<- "scrsurv"
  return(sout)
  
}


                      
scrsurvfit<- function(time1,event1,time2,event2,copulaparam,copulafam,
                      copulaparams = NULL,bootind =NULL,zetas =NULL,
                      model,casewt,surv2km,method,confint, se.method,B,
                      conftype){
  zstar<- qnorm(0.05/2,lower.tail = FALSE)
  n<- length(time1)
  if(method=="J"){
    sout<- scrsurvfit.JFKC(time1 = time1,event1 = event1,time2 = time2,event2 = event2,
                           copulaparam = copulaparam,copulafam = copulafam,
                           model = model,casewt=casewt,surv2km=surv2km)
    
    names(sout)<- c("t1.surv","t2.surv")
    if(isTRUE(confint)){
      S1list<- list()
      H1list<- list()
      S2list<- list()
      H2list<- list()
      cat("Please wait for a while...")
      if(conftype %in%1:2){
        if(se.method =="resampling"){
          
          for(k in 1:B){
            zeta<-rexp(n,rate=1)
            fit<- scrsurvfit.JFKC(time1 = time1,event1 = event1,time2 = time2,event2 = event2,
                                  copulaparam = copulaparam,copulafam = copulafam,
                                  model = model,casewt=casewt*zeta,surv2km=surv2km)
            
            S1list[[k]]<- stepfun(fit$s1$time,c(1,fit$s1$surv),right = FALSE)(sout$t1.surv$time)
            H1list[[k]]<- stepfun(fit$s1$time,c(0,fit$s1$cumhaz),right = FALSE)(sout$t1.surv$time)
            if(!isTRUE(surv2km)){
              S2list[[k]]<- stepfun(fit$s2$time,c(1,fit$s2$surv),right = FALSE)(sout$t2.surv$time)
              H2list[[k]]<- stepfun(fit$s2$time,c(0,fit$s2$cumhaz),right = FALSE)(sout$t2.surv$time)
            }
            if(k%%20==0){cat(round(k/B,1)*100,"%.")}
          }
          
        }else{
          for(k in 1:B){
            indd<- sort(sample(1:n,replace = TRUE))
            fit<- scrsurvfit.JFKC(time1 = time1[indd],event1 = event1[indd],
                                  time2 = time2[indd],event2 = event2[indd],
                                  copulaparam = copulaparam,copulafam = copulafam,
                                  model = model,casewt=casewt[indd],surv2km=surv2km)
            
            S1list[[k]]<- stepfun(fit$s1$time,c(1,fit$s1$surv),right = FALSE)(sout$t1.surv$time)
            H1list[[k]]<- stepfun(fit$s1$time,c(0,fit$s1$cumhaz),right = FALSE)(sout$t1.surv$time)
            if(!isTRUE(surv2km)){
              S2list[[k]]<- stepfun(fit$s2$time,c(1,fit$s2$surv),right = FALSE)(sout$t2.surv$time)
              H2list[[k]]<- stepfun(fit$s2$time,c(0,fit$s2$cumhaz),right = FALSE)(sout$t2.surv$time)
            }
            if(k%%20==0){cat(round(k/B,1)*100,"%.")}
          }
        }
        
        S1list<- do.call(rbind,S1list)
        H1list<- do.call(rbind,H1list)
        sout$t1.surv$std.err<- apply(S1list,2,sd,na.rm=TRUE)
        sout$t1.surv$std.chaz<- apply(H1list,2,sd,na.rm=TRUE)
        if(conftype==2){
          sout$t1.surv$lower<- apply(S1list,2,quantile,probs = 0.05/2,na.rm=TRUE)
          sout$t1.surv$upper<- apply(S1list,2,quantile,probs = 1-0.05/2,na.rm=TRUE)
        }else if(conftype==1){
          sout$t1.surv$lower<- pmin(1,pmax(apply(S1list,2,mean,na.rm=TRUE)-
                                             zstar*sout$t1.surv$std.err,0))
          sout$t1.surv$upper<- pmin(1,pmax(apply(S1list,2,mean,na.rm=TRUE)+
                                             zstar*sout$t1.surv$std.err,0))
        }
        
        if(!isTRUE(surv2km)){
          S2list<- do.call(rbind,S2list)
          H2list<- do.call(rbind,H2list)
          sout$t2.surv$std.err<- apply(S2list,2,sd,na.rm=TRUE)
          sout$t2.surv$std.chaz<- apply(H2list,2,sd,na.rm=TRUE)
          if(conftype==2){
            sout$t2.surv$lower<- apply(S2list,2,quantile,probs = 0.05/2,na.rm=TRUE)
            sout$t2.surv$upper<- apply(S2list,2,quantile,probs = 1-0.05/2,na.rm=TRUE)
          }else if(conftype==1){
            sout$t2.surv$lower<- pmin(1,pmax(apply(S2list,2,mean,na.rm=TRUE)-
                                               zstar*sout$t2.surv$std.err,0))
            sout$t2.surv$upper<- pmin(1,pmax(apply(S2list,2,mean,na.rm=TRUE)+
                                               zstar*sout$t2.surv$std.err,0))
          }
          
        }
      }
      if(conftype %in%3:4){
        if(se.method =="resampling"){
          
          for(k in 1:B){
            zeta<- zetas[k,]
            fit<- scrsurvfit.JFKC(time1 = time1,event1 = event1,time2 = time2,event2 = event2,
                                  copulaparam = copulaparams[k],copulafam = copulafam,
                                  model = model,casewt=casewt*zeta,surv2km=surv2km)
            
            S1list[[k]]<- stepfun(fit$s1$time,c(1,fit$s1$surv),right = FALSE)(sout$t1.surv$time)
            H1list[[k]]<- stepfun(fit$s1$time,c(0,fit$s1$cumhaz),right = FALSE)(sout$t1.surv$time)
            if(!isTRUE(surv2km)){
              S2list[[k]]<- stepfun(fit$s2$time,c(1,fit$s2$surv),right = FALSE)(sout$t2.surv$time)
              H2list[[k]]<- stepfun(fit$s2$time,c(0,fit$s2$cumhaz),right = FALSE)(sout$t2.surv$time)
            }
            if(k%%20==0){cat(round(k/B,1)*100,"%.")}
          }
          
        }else{
          for(k in 1:B){
            indd<- bootind[k,]
            fit<- scrsurvfit.JFKC(time1 = time1[indd],event1 = event1[indd],
                                  time2 = time2[indd],event2 = event2[indd],
                                  copulaparam = copulaparams[k],copulafam = copulafam,
                                  model = model,casewt=casewt[indd],surv2km=surv2km)
            
            S1list[[k]]<- stepfun(fit$s1$time,c(1,fit$s1$surv),right = FALSE)(sout$t1.surv$time)
            H1list[[k]]<- stepfun(fit$s1$time,c(0,fit$s1$cumhaz),right = FALSE)(sout$t1.surv$time)
            if(!isTRUE(surv2km)){
              S2list[[k]]<- stepfun(fit$s2$time,c(1,fit$s2$surv),right = FALSE)(sout$t2.surv$time)
              H2list[[k]]<- stepfun(fit$s2$time,c(0,fit$s2$cumhaz),right = FALSE)(sout$t2.surv$time)
            }
            if(k%%20==0){cat(round(k/B,1)*100,"%.")}
          }
        }
        
        S1list<- do.call(rbind,S1list)
        H1list<- do.call(rbind,H1list)
        sout$t1.surv$std.err<- apply(S1list,2,sd,na.rm=TRUE)
        sout$t1.surv$std.chaz<- apply(H1list,2,sd,na.rm=TRUE)
        if(conftype==4){
          sout$t1.surv$lower<- apply(S1list,2,quantile,probs = 0.05/2,na.rm=TRUE)
          sout$t1.surv$upper<- apply(S1list,2,quantile,probs = 1-0.05/2,na.rm=TRUE)
        }else if(conftype==3){
          sout$t1.surv$lower<- pmin(1,pmax(apply(S1list,2,mean,na.rm=TRUE)-
                                             zstar*sout$t1.surv$std.err,0))
          sout$t1.surv$upper<- pmin(1,pmax(apply(S1list,2,mean,na.rm=TRUE)+
                                             zstar*sout$t1.surv$std.err,0))
        }
        
        if(!isTRUE(surv2km)){
          S2list<- do.call(rbind,S2list)
          H2list<- do.call(rbind,H2list)
          sout$t2.surv$std.err<- apply(S2list,2,sd,na.rm=TRUE)
          sout$t2.surv$std.chaz<- apply(H2list,2,sd,na.rm=TRUE)
          if(conftype==4){
            sout$t2.surv$lower<- apply(S2list,2,quantile,probs = 0.05/2,na.rm=TRUE)
            sout$t2.surv$upper<- apply(S2list,2,quantile,probs = 1-0.05/2,na.rm=TRUE)
          }else if(conftype==3){
            sout$t2.surv$lower<- pmin(1,pmax(apply(S2list,2,mean,na.rm=TRUE)-
                                               zstar*sout$t2.surv$std.err,0))
            sout$t2.surv$upper<- pmin(1,pmax(apply(S2list,2,mean,na.rm=TRUE)+
                                               zstar*sout$t2.surv$std.err,0))
          }
          
        }
      }
    }
    
  }
  
  
  if(method!="J"){
    S2fitkm<- survival::survfit(Surv(time2,event2)~1,model = model) 
    
    
    S1fit<- switch (method,
                    "L" = scrsurvfit.LRA(time1 = time1,event1 = event1,time2 = time2,event2 = event2,
                                         copulaparam = copulaparam,copulafam = copulafam,
                                         model = model,casewt=casewt),
                    "F" = scrsurvfit.FJC(time1 = time1,event1 = event1,time2 = time2,event2 = event2,
                                         copulaparam = copulaparam,copulafam = copulafam,
                                         model = model,casewt=casewt,S2=S2fitkm),
                    "K" = scrsurvfit.kernel(time1 = time1,event1 = event1,time2 = time2,event2 = event2,S2=S2fitkm)
    )
    
    
    if(isTRUE(confint)){
      S1list<- list()
      H1list<- list()
      
      if(conftype %in%1:2){
        if(method=="L"){
          if(se.method=="jackknife"){
            for (k in 1:n) {
              fit<- scrsurvfit.LRA(time1 = time1[-k],event1 = event1[-k],time2 = time2[-k],event2 = event2[-k],
                                   copulaparam = copulaparam,copulafam = copulafam,model = model,casewt=casewt[-k])
              S1list[[k]]<- stepfun(fit$time,c(1,fit$surv),right = FALSE)(S1fit$time)
              H1list[[k]]<- stepfun(fit$time,c(0,fit$cumhaz),right = FALSE)(S1fit$time)
            }
          }else{
            for (k in 1:B) {
              indd<- sort(sample(1:n,replace = TRUE))
              fit<- scrsurvfit.LRA(time1 = time1[indd],event1 = event1[indd],
                                   time2 = time2[indd],event2 = event2[indd],
                                   copulaparam = copulaparam,copulafam = copulafam,
                                   model = model,casewt=casewt[indd])
              S1list[[k]]<- stepfun(fit$time,c(1,fit$surv),right = FALSE)(S1fit$time)
              H1list[[k]]<- stepfun(fit$time,c(0,fit$cumhaz),right = FALSE)(S1fit$time)
            }
          }
          
          S1list<- do.call(rbind,S1list)
          H1list<- do.call(rbind,H1list)
          S1fit$std.err<- apply(S1list,2,sd,na.rm=TRUE)
          S1fit$std.chaz<- apply(H1list,2,sd,na.rm=TRUE)
          if(conftype==2){
            S1fit$lower<- apply(S1list,2,quantile,probs = 0.05/2,na.rm=TRUE)
            S1fit$upper<- apply(S1list,2,quantile,probs = 1-0.05/2,na.rm=TRUE)
          }else if(conftype==1){
            S1fit$lower<- pmin(1,pmax(apply(S1list,2,mean,na.rm=TRUE)-
                                        zstar*S1fit$std.err,0))
            S1fit$upper<- pmin(1,pmax(apply(S1list,2,mean,na.rm=TRUE)+
                                        zstar*S1fit$std.err,0))
          }
          
        }else{
          
          if(se.method=="jackknife"){
            for (k in 1:n) {
              
              fit2<- survival::survfit(Surv(time2[-k],event2[-k])~1,model = model) 
              
              fit1<- switch (method,
                             "F" = scrsurvfit.FJC(time1 = time1[-k],event1 = event1[-k],time2 = time2[-k],event2 = event2[-k],
                                                  copulaparam = copulaparam,copulafam = copulafam,
                                                  model = model,casewt=casewt[-k],S2=fit2),
                             "K" = scrsurvfit.kernel(time1 = time1[-k],event1 = event1[-k],time2 = time2[-k],
                                                     event2 = event2[-k], S2=fit2)
              )
              S1list[[k]]<- stepfun(fit1$time,c(1,fit1$surv),right = FALSE)(S1fit$time)
              H1list[[k]]<- stepfun(fit1$time,c(0,fit1$cumhaz),right = FALSE)(S1fit$time)
              
            }
          }else{
            for (k in 1:B) {
              indd<- sort(sample(1:n,replace = TRUE))
              fit2<- survival::survfit(Surv(time2[indd],event2[indd])~1,model = model) 
              
              fit1<- switch (method,
                             "F" = scrsurvfit.FJC(time1 = time1[indd],event1 = event1[indd],time2 = time2[indd],event2 = event2[indd],
                                                  copulaparam = copulaparam,copulafam = copulafam,
                                                  model = model,casewt=casewt[indd],S2=fit2),
                             "K" = scrsurvfit.kernel(time1 = time1[indd],event1 = event1[indd],time2 = time2[indd],
                                                     event2 = event2[indd],S2=fit2)
              )
              S1list[[k]]<- stepfun(fit1$time,c(1,fit1$surv),right = FALSE)(S1fit$time)
              H1list[[k]]<- stepfun(fit1$time,c(0,fit1$cumhaz),right = FALSE)(S1fit$time)
              
            }
          }
          
          S1list<- do.call(rbind,S1list)
          H1list<- do.call(rbind,H1list)
          S1fit$std.err<- apply(S1list,2,sd,na.rm=TRUE)
          S1fit$std.chaz<- apply(H1list,2,sd,na.rm=TRUE)
          if(conftype==2){
            S1fit$lower<- apply(S1list,2,quantile,probs = 0.05/2,na.rm=TRUE)
            S1fit$upper<- apply(S1list,2,quantile,probs = 1-0.05/2,na.rm=TRUE)
          }else if(conftype==1){
            S1fit$lower<- pmin(1,pmax(apply(S1list,2,mean,na.rm=TRUE)-
                                        zstar*S1fit$std.err,0))
            S1fit$upper<- pmin(1,pmax(apply(S1list,2,mean,na.rm=TRUE)+
                                        zstar*S1fit$std.err,0))
          }
          
        }
      }
      if(conftype %in%3:4){
        if(method=="L"){
          for (k in 1:B) {
            indd<-  bootind[k,]
            fit<- scrsurvfit.LRA(time1 = time1[indd],event1 = event1[indd],
                                 time2 = time2[indd],event2 = event2[indd],
                                 copulaparam = copulaparams[k],copulafam = copulafam,
                                 model = model,casewt=casewt[indd])
            S1list[[k]]<- stepfun(fit$time,c(1,fit$surv),right = FALSE)(S1fit$time)
            H1list[[k]]<- stepfun(fit$time,c(0,fit$cumhaz),right = FALSE)(S1fit$time)
          }
          
          S1list<- do.call(rbind,S1list)
          H1list<- do.call(rbind,H1list)
          S1fit$std.err<- apply(S1list,2,sd,na.rm=TRUE)
          S1fit$std.chaz<- apply(H1list,2,sd,na.rm=TRUE)
          if(conftype==4){
            S1fit$lower<- apply(S1list,2,quantile,probs = 0.05/2,na.rm=TRUE)
            S1fit$upper<- apply(S1list,2,quantile,probs = 1-0.05/2,na.rm=TRUE)
          }else if(conftype==3){
            S1fit$lower<- pmin(1,pmax(apply(S1list,2,mean,na.rm=TRUE)-
                                        zstar*S1fit$std.err,0))
            S1fit$upper<- pmin(1,pmax(apply(S1list,2,mean,na.rm=TRUE)+
                                        zstar*S1fit$std.err,0))
          }
          
        }else{
          
          for (k in 1:B) {
            indd<-  bootind[k,]
            fit2<- survival::survfit(Surv(time2[indd],event2[indd])~1,model = model) 
            
            fit1<- switch (method,
                           "F" = scrsurvfit.FJC(time1 = time1[indd],event1 = event1[indd],time2 = time2[indd],event2 = event2[indd],
                                                copulaparam = copulaparams[k],copulafam = copulafam,
                                                model = model,casewt=casewt[indd],S2=fit2),
                           "K" = scrsurvfit.kernel(time1 = time1[indd],event1 = event1[indd],time2 = time2[indd],
                                                   event2 = event2[indd],S2=fit2)
            )
            S1list[[k]]<- stepfun(fit1$time,c(1,fit1$surv),right = FALSE)(S1fit$time)
            H1list[[k]]<- stepfun(fit1$time,c(0,fit1$cumhaz),right = FALSE)(S1fit$time)
            
          }
          
          
          S1list<- do.call(rbind,S1list)
          H1list<- do.call(rbind,H1list)
          S1fit$std.err<- apply(S1list,2,sd,na.rm=TRUE)
          S1fit$std.chaz<- apply(H1list,2,sd,na.rm=TRUE)
          if(conftype==4){
            S1fit$lower<- apply(S1list,2,quantile,probs = 0.05/2,na.rm=TRUE)
            S1fit$upper<- apply(S1list,2,quantile,probs = 1-0.05/2,na.rm=TRUE)
          }else if(conftype==3){
            S1fit$lower<- pmin(1,pmax(apply(S1list,2,mean,na.rm=TRUE)-
                                        zstar*S1fit$std.err,0))
            S1fit$upper<- pmin(1,pmax(apply(S1list,2,mean,na.rm=TRUE)+
                                        zstar*S1fit$std.err,0))
          }
          
        }
      }
      
    }
    
    sout<- list(t1.surv=S1fit,t2.surv=S2fitkm)
  }
  return(sout)
}




scrsurvfit.FJC<- function(time1,event1,time2,event2,copulaparam,copulafam,model,casewt,S2){
  
  ptime<- unique(as.numeric(sort(time1)))
  ptime2<- unique(as.numeric(sort(time1[event1==1])))
  
  t2surv<- approx(c(0,S2$time),c(1,S2$surv),yright = min(S2$surv),ties = mean,
                  method = "constant",xout = ptime2,f = 0)$y
  
  
  event3<- event1+(1-event1)*event2
  t12fit <- prodlim::prodlim(Surv(time1,event3)~1)
  
  t12surv<- predict(t12fit,times = ptime2,level.chaos='1',type="surv",mode="matrix",bytime=TRUE)
  
  
  surv1<- sapply(archmCopulaLink(copulafam,param = copulaparam, p=t12surv)-
                   archmCopulaLink(copulafam,param = copulaparam, p=t2surv),
                 FUN = archmCopulaLink_inv,copulafam = copulafam,param = copulaparam)
  surv1[surv1>1]<- 1
  surv1[surv1<0]<- 0
  
  F1<- quantreg::rearrange(stepfun(x= ptime2,y = c(0,1-surv1),right = FALSE,ties = mean))(ptime)
  surv1<- 1-F1 
  
  cumhaz<-  -log(surv1)
  cumhaz[cumhaz<0 ]<- 0
  
  
  n.risk<- as.numeric(sapply(ptime,function(tt)sum(time1>=tt)))
  n.event<- as.numeric(sapply(ptime,function(tt)sum(event1*(time1==tt))))
  n.censor<- as.numeric(sapply(ptime,function(tt)sum((1-event1)*(time1==tt))))
  
  
  out<- list(time=ptime,n.risk=n.risk,n.event=n.event,n.censor=n.censor,
               surv=surv1,cumhaz=cumhaz,type="right")
  
  
  return(out)
  
}



scrsurvfit.LRA<- function(time1,event1,time2,event2,copulaparam,copulafam,
                         model,casewt){
  
  ptime<- as.numeric(sort(time1))
  ptime2<- unique(as.numeric(sort(time1[event1==1])))
  
  event3<- event1+(1-event1)*event2
  t12fit <- prodlim::prodlim(Surv(time1,event3)~1)
  
  t12surv<- predict(t12fit,times = ptime,level.chaos='1',type="surv",mode="matrix",bytime=TRUE)
  
  surv1<- sapply(ptime2,function(t){
    pos<- max(which(ptime==t))
    archmCopulaLink_inv(copulafam,param = copulaparam,
                        y=sum(archmCopulaLink(copulafam,param = copulaparam, p=t12surv[1:pos])-
                                archmCopulaLink(copulafam,param = copulaparam, p=c(1,t12surv)[1:pos]),
                              na.rm = TRUE))
  })
  
  surv1[surv1>1]<- 1
  surv1[surv1<0]<- 0
  
  ptime<- unique(ptime)
  
  F1<- quantreg::rearrange(stepfun(x= ptime2,y = c(0,1-surv1),right = FALSE,ties = mean))(ptime)
  surv1<- 1-F1 
  
  cumhaz<-  -log(surv1)
  cumhaz[cumhaz<0 ]<- 0
  
  
  
  n.risk<- as.numeric(sapply(ptime,function(tt)sum(time1>=tt)))
  n.event<- as.numeric(sapply(ptime,function(tt)sum(event1*(time1==tt))))
  n.censor<- as.numeric(sapply(ptime,function(tt)sum((1-event1)*(time1==tt))))
  
  
  out<- list(time=ptime,n.risk=n.risk,n.event=n.event,n.censor=n.censor,
               surv=surv1,cumhaz=cumhaz,type="right")
  
  
  
  return(out)
  
}

scrsurvfit.JFKC<-function(time1,event1,time2,event2,copulaparam,copulafam,
                          model,casewt,surv2km=FALSE){
  
  tp<- sort(unique(c(time1[event1==1],time2[event2==1])))
  # tp<- sort(unique(c(time1,time2)))
  term1 <- sapply(tp,function(t) sum((time1> t)*casewt,na.rm=TRUE))/sum(casewt)
  if(!isTRUE(surv2km)){
    term2 <- sapply(tp,function(t) sum((time2> t)*casewt,na.rm=TRUE))/sum(casewt)
  }
  
  S1<- survival::survfit(Surv(time1,event1)~1,weights = casewt)
  F1bar.ini<- stepfun(x = c(0,S1$time),y= c(1,S1$surv,0),right = FALSE,ties = mean)(tp)
  
  S2<- survival::survfit(Surv(time2,event2)~1,weights = casewt)
  F2bar.ini<-  stepfun(x = c(0,S2$time),y= c(1,S2$surv,0),right = FALSE,ties = mean)(tp)
  
  it<- 1
  mit<- 200
  test0<- 1000
  convergence<- 0
  repeat{
    F1bar.ini.x<- approx(x=tp,y=F1bar.ini,xout = time1,method = "constant",f = 0,ties = mean)$y
    F2bar.ini.x<- approx(x=tp,y=F2bar.ini,xout = time2,method = "constant",f = 0,ties = mean)$y
    
    dc10<- sapply(1:length(tp),function(k)
      sum(Copulafn(copulafam = copulafam,param = copulaparam,p1 = F1bar.ini[k],p2 = F2bar.ini.x)/
            Copulafn(copulafam = copulafam,param = copulaparam,p1 = F1bar.ini.x,p2 = F2bar.ini.x)*
            (1-event1)*(1-event2)*casewt*(time1<= tp[k]),na.rm = TRUE)/sum(casewt))
    
    dc12<- sapply(1:length(tp),function(k)
      sum(dev_Copula(copulafam = copulafam,param = copulaparam,p1 = F1bar.ini[k],p2 = F2bar.ini.x,mode = "2")/
            dev_Copula(copulafam = copulafam,param = copulaparam,p1 = F1bar.ini.x,p2 = F2bar.ini.x,mode = "2")*
            (1-event1)*event2*casewt*(time1<= tp[k]),na.rm = TRUE)/sum(casewt))
    dc10[is.na(dc10)]<- 0
    dc10[is.infinite(dc10)]<- 0
    
    dc12[is.na(dc12)]<- 0
    dc12[is.infinite(dc12)]<- 0
    F1bar.est<- term1 + dc10+ dc12
    
    if(!isTRUE(surv2km)){
      dc20<- sapply(1:length(tp),function(k)
        sum(Copulafn(copulafam = copulafam,param = copulaparam,p1 = F1bar.ini.x,p2 = F2bar.ini[k])/
              Copulafn(copulafam = copulafam,param = copulaparam, p1 = F1bar.ini.x,p2 = F2bar.ini.x)*
              (1-event1)*(1-event2)*casewt*(time2<= tp[k]),na.rm = TRUE)/sum(casewt))
      
      dc21<- sapply(1:length(tp),function(k)
        sum(dev_Copula(copulafam = copulafam,param = copulaparam,p1 = F1bar.ini.x,p2 = F2bar.ini[k],mode = "1")/
              dev_Copula(copulafam = copulafam,param = copulaparam, p1 = F1bar.ini.x,p2 = F2bar.ini.x,mode = "1")*
              event1*(1-event2)*casewt*(time2<= tp[k]),na.rm = TRUE)/sum(casewt))
      dc20[is.na(dc20)]<- 0
      dc20[is.infinite(dc20)]<- 0
      
      dc21[is.na(dc21)]<- 0
      dc21[is.infinite(dc21)]<- 0
      F2bar.est<- term2 + dc20 + dc21
      F2bar.est<- pmin(pmax(0,F2bar.est),1)
    }else{
      F2bar.est<-F2bar.ini
    }
    
    F1bar.est<- pmin(pmax(0,F1bar.est),1)
    
    
    test<- max(max(abs(F1bar.est-F1bar.ini)),max(abs(F2bar.est-F2bar.ini)))
    # if(it%%20==1){cat(it,test,".\n")}
    if(test< 1e-4|it>mit){
      # cat(it,test,".\n")
      convergence=0
      break
    }else if(test>test0){
      # cat(it,test,".\n")
      F1bar.est <- F1bar.ini
      F2bar.est <- F2bar.ini
      if(test>0.05){convergence=1}
      break
    }else{
      it<-it+1
      F1bar.ini<- F1bar.est
      F2bar.ini<- F2bar.est
      test0<- test
    }
    
  }
  
  
  
  F1<- 1- F1bar.est
  F2<- 1- F2bar.est
  
  napos<- union(which(is.na(F1)),which(is.na(F2)))
  if(length(napos)>0){
    tp<- tp[-napos]
    F1<- F1[-napos]
    F2<- F2[-napos]
    
  }
  
  ptime1<-  as.numeric(sort(unique(time1)))
  ptime2<-  as.numeric(sort(unique(time2)))
  F1<- quantreg::rearrange(stepfun(x= tp,y = c(0,F1),right = FALSE,ties = mean))(ptime1)
  F2<- quantreg::rearrange(stepfun(x= tp,y = c(0,F2),right = FALSE,ties = mean))(ptime2)
  
  
  surv1<- 1-F1 
  surv2<- 1-F2
  
  n.risk1<- as.numeric(sapply(ptime1,function(tt)sum(time1>=tt)))
  n.event1<- as.numeric(sapply(ptime1,function(tt)sum(event1*(time1==tt))))
  n.censor1<- as.numeric(sapply(ptime1,function(tt)sum((1-event1)*(time1==tt))))
  
  cumhaz1<-  -log(surv1)
  cumhaz1[cumhaz1<0 ]<- 0
  
  
  
  out1<- list(time=ptime1,n.risk=n.risk1,n.event=n.event1,n.censor=n.censor1,
              surv=surv1,cumhaz=cumhaz1,type="right")
  if(!isTRUE(surv2km)){
    n.risk2<- as.numeric(sapply(ptime2,function(tt)sum(time2>=tt)))
    n.event2<- as.numeric(sapply(ptime2,function(tt)sum(event2*(time2==tt))))
    n.censor2<- as.numeric(sapply(ptime2,function(tt)sum((1-event2)*(time2==tt))))
    
    cumhaz2<-  -log(surv2)
    cumhaz2[cumhaz2<0 ]<- 0
    out2<- list(time=ptime2,n.risk=n.risk2,n.event=n.event2,n.censor=n.censor2,
                surv=surv2,cumhaz=cumhaz2,type="right")
  }else{
    out2<- S2
  }
  
  return(list(s1=out1,s2=out2))
  
  
  
  
}
scrsurvfit.kernel<- function(time1,event1,time2,event2,S2){
  # A modification of codes from Nevo D, Gorfine M. (2022)
  # Causal inference for semi-competing risks data. 
  # Biostatistics, 14, 23(4):1115-1132.
  ptime<- unique(as.numeric(sort(time1)))
  n.times<- length(ptime)
  
  dS2 <-  diff(c(1, S2$surv))
  dS2.times <-  S2$time
  dS2.times<- dS2.times[dS2!= 0]
  dS2 <- dS2[dS2!= 0]
  n.times2 <- length(dS2.times)
  
  T1dead <- time1[event2==1]
  event1dead <- event1[event2==1]
  T2dead <- time2[event2==1]
  
  # Calculate S_{1|A=0, T_2=t}
  S1_2 <- prodlim::prodlim(Surv(T1dead, event1dead) ~ T2dead)
  newdata <- data.frame(T2dead = dS2.times)
  S1_2 <- predict(S1_2, newdata = newdata, times =  dS2.times)
  S1_2 <- do.call(rbind, S1_2)
  
  
  # Change NA after no more events into the last estimated value
  S1_2 <- t(apply(S1_2, 1, RepNAmin))
  
  S1.ptime <-  vector(length = n.times)
  S1st <- matrix(nrow = n.times, ncol = n.times2)
  
  S1AtT2s <- function(my.t, my.s) {
    if (my.s < my.t){      my.t <- my.s}
    s.place <- findInterval(my.s, dS2.times)
    t.place <- findInterval(my.t, dS2.times)
    if (t.place==0 | s.place==0) {
      S.res <- 1
    } else {
      t1.curve.s <- S1_2[s.place, ]
      if (t.place==length(dS2.times)) {
        S.res <-  t1.curve.s[length(dS2.times)]
      } else {
        S.res <- t1.curve.s[t.place]
      }}
    
    return(S.res)
  }
  for (j in 1:n.times) {
    t.now <- ptime[j]
    pos.t.now <- findInterval(t.now, dS2.times)
    
    if(pos.t.now==0) {
      S1.ptime[j] <- 1
    } else {
      S1st[j, ] <-  sapply(dS2.times, S1AtT2s, my.t = t.now)
      S1.ptime[j] <- -sum(S1st[j, ] * dS2)
    }
    
    
    
  }
  
  
  surv1<- S1.ptime
  surv1[surv1>1]<- 1
  surv1[surv1<0]<- 0
  F1<- quantreg::rearrange(stepfun(x= ptime,y = c(0,1-surv1),right = FALSE,ties = mean))(ptime)
  surv1<- 1-F1 
  
  cumhaz<-  -log(surv1)
  cumhaz[cumhaz<0 ]<- 0
  
  
  n.risk<- as.numeric(sapply(ptime,function(tt)sum(time1>=tt)))
  n.event<- as.numeric(sapply(ptime,function(tt)sum(event1*(time1==tt))))
  n.censor<- as.numeric(sapply(ptime,function(tt)sum((1-event1)*(time1==tt))))
  
  
  out<- list(time=ptime,n.risk=n.risk,n.event=n.event,n.censor=n.censor,
             surv=surv1,cumhaz=cumhaz,type="right")
  
  return(out)
}

RepNAmin <- function(x){
  min.x <- min(x, na.rm = T)
  x[is.na(x)] <- min.x
  return(x)
}

