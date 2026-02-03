
#' @title  Simulating a semi-competing risks data with multiple intermediate events
#' @aliases simSCRmul
#' @description
#' Simulating a semi-competing risks data with multiple intermediate events from a copula-based parametric model
#' with the Archimedean copula structure.
#'
#' @usage simSCRmul(n  = 100, cens.rate = 20, K = 7, copulafam = "frank",
#'   tau.alpha = 0.5, tau.theta = seq(0.85, 0.2,length.out = K),
#'   params = list(rate1 = 1, rate2 = 0.6))
#'
#' @param n Sample size
#' @param cens.rate Parameter for the censorship distribution.
#' The censoring variable (C) is generated from a \eqn{Uniform[0,cens.rate]} distribution
#' @param K The number of intermediate events. \code{K} is an integer and greater than 1.
#' @param copulafam a character string specifying the family of an Archimedean copula.
#' Currently supported families are "frank", "clayton",  "amh", "gumbel", and "joe".
#' @param tau.alpha Kendall's tau corresponding to copula parameter \eqn{\alpha}.
#' @param tau.theta Kendall's tau corresponding to copula parameter \eqn{\theta}.
#' (See more in Details)
#' @param params a list of parameters rate1 and rate2.
#' The marginal distributions of intermediate events are identically \eqn{Exp(rate1)},
#' and the marginal of the terminal event time is \eqn{Exp(rate2)}.
#'
#' @return A
#' @export
#' @seealso \code{\link{mscr}}
#'
simSCRmul<- function (n  = 100,cens.rate = 20,K=7,
                      copulafam="frank",tau.alpha=0.5,
                      tau.theta=seq(0.85,0.2,length.out=K),
                     params=list(rate1=1, rate2=0.6))
{
  call<- match.call()
  if(length(tau.theta)==1){tau.theta<- rep(tau.theta,K)}
  alpha<- Calitau(copulafam,tau.alpha)
  theta<- sapply(tau.theta,Calitau,copulafam=copulafam)
  death<- stats::rexp(n =n,rate = params$rate2)
  surv.death<- stats::pexp(death,rate = params$rate2,lower.tail = FALSE)

  copula_type <- copula::archmCopula(family= copulafam,param = alpha,dim = K)
  U <- copula::rCopula(n, copula_type) #Generate uniform random numbers from the copula

  # Transform uniform random numbers to samples from the desired distribution
  # Apply inverse transform method using marginal CDFs
  intermediate <- u.inter<- matrix(NA, nrow = n, ncol = K)
  for (k in 1:K) {
    for (i in 1:n) {
      u.inter[i,k]<- stats::uniroot(f = function(xx){
        mycsurv(p1=xx,p2=surv.death[i],copulafam=copulafam,copulaparam=theta[k])-
          U[i, k]}, interval = c(0, 1))$root
    }



    intermediate[, k] <- stats::qexp(u.inter[,k],rate = params$rate1,lower.tail = FALSE) #  inverse of marginal CDF
  }
  colnames(intermediate)<- paste0("T",1:k)

  C <- stats::runif(n,min = 0 ,max = cens.rate)

  T2obs<- pmin(death,C)
  T1obs<- apply(intermediate,2,function(xx)pmin(xx,T2obs))
  delta<- 1*(T1obs==intermediate)
  delta.tilde<- 1*(death< C)

  data<- list(T1= T1obs,T2= T2obs,event1= delta, event2= delta.tilde)
  # class(data)<- "mscrData"
  return(data)

}

## #' Return the First Parts of the mscrData object.
## #' @usage head(x, ...)
## #' @param x the mscrData object
## #' @return A list of 2 length with components representing observations of intermediate event times
## #' and terminal event time, respectively.
## #' @export

# "head.mscrData"<- function(x,...){
#
#   h<- list()
#   h[[1]]<- head(cbind(mdata$T1,mdata$event1),...)
#   h[[2]]<- head(cbind(mdata$T2,mdata$event2),...)
#   colnames(h[[2]])<- c("DTH","status")
#   names(h)<-  c("intermediate event times","terminal event time")
#
#   h
#
# }





#' @noRd
mycsurv<- function(p1,p2,copulafam,copulaparam){
  out<- dev_Copula(copulafam = copulafam,param = copulaparam,p1 =p1,p2 =p2,mode = "2")
  return(out)

}
