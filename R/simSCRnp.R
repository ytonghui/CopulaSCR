###### simulate a random sample ########
#' @title Simulating a semi-competing risks data from a copula-based parametric model
#' @aliases simSCR
#' @usage simSCR(n  = 100, cens.rate = 5, copulafam = "frank", tau = 0.5,
#' params = list(marginsDist = rep("exp", 2),rate1 = 1,rate2 = 1))
#' @description
#' Simulating a semi-competing risks data from a copula-based parametric model
#' with the Archimedean copula structure.
#'
#' @param n Sample size
#' @param cens.rate Parameter for the censorship distribution.
#' The censoring variable (C) is generated from a \eqn{Uniform[0,cens.rate]} distribution.
#' @param copulafam a character string specifying the family of an Archimedean copula.
#' Currently supported families are "frank", "clayton",  "amh", "gumbel", and "joe".
#' @param tau  Kendall's tau of an Archimedean copula.
#' @param params a list of marginsDist, rate1, and rate2. marginsDist is object of class "character" specifying the marginal distributions of nonterminal (T1) and terminal (T2) event times in sequence.
#' rate1 and rate2 are parameter values of their marginal distributions.
#'
#' @return A data.frame containing semi-competing risks outcomes from \code{n} subjects.
#' It is of dimension \eqn{n\times 4} with the columns named as T1, T2, event1 and event2. \cr
#' \describe{
#'   \item{T1}{a vector of the observed nonterminal event times for the n individuals, i.e.,
#'    \eqn{\min(\tilde{T}_1,\tilde{T}_2,C)}, where \eqn{\tilde{T}_1}, and \eqn{\tilde{T}_2}
#'    are the true nonterminal and terminal event times, respectively.}
#'   \item{T2}{a vector of the observed terminal event times for the n individuals, i.e.,
#'    \eqn{\tilde{T}_2,C)}.}
#'   \item{event1}{a vector of censoring indicators for the non-terminal event time (1=event occurred, 0=censored).}
#'   \item{event2}{a vector of censoring indicators for the terminal event time (1=event occurred, 0=censored).}
#' }
## #' @import copula
#'
#' @export
#' @seealso \code{\link{scrassonp}}, \code{\link{scrsurv}}



simSCR<- function (n  = 100,cens.rate = 5,copulafam="frank",tau=0.5,
                     params=list(marginsDist = rep("exp",2),rate1=1,rate2=1))
{
  alpha<-  Calitau(copulafam = copulafam,tau = tau)

  biSurvTimes<-  copula::rMvdc(n,copula::mvdc(copula =copula::archmCopula(family= copulafam,
                                                  param = alpha,
                                                  dim = 2),
                              margins = params$marginsDist,
                              paramMargins = c(list(list(rate = params$rate1)),
                                               list(list(rate = params$rate2)))))
  T1<- biSurvTimes[,1]
  T2<- biSurvTimes[,2]

  C <- stats::runif(n,min = 0 ,max = cens.rate)

  T2obs<- pmin(T2,C)
  T1obs<- pmin(T1,T2obs)
  delta1<- 1*(T1<= T2obs)
  delta2<- 1*(T2<= C)

  cat("observed rate:", sum(delta1)/n,"for T1,",sum(delta2)/n, "for T2","\n")


  data<- data.frame(T1= T1obs,T2= T2obs,event1= delta1, event2= delta2)
  return(data)

}



#' @title Simulating a semi-competing risks data from a copula-based parametric model
#' @aliases simSCRtr
#' @usage simSCRtr(n  = 100, K = 2, cens.rate = 5, copulafam = "frank", tau = 0.5,
#' params = list(marginsDist = rep("exp", 2),rate1 = 1,rate2 = 1))
#'
#' @description
#' Simulating a stratified semi-competing risks data from a copula-based parametric model
#' with the Archimedean copula structure.
#'
#' @param n Sample size
#' @param K the number of types of treatments
#' @param cens.rate Parameter for the censorship distribution.
#' The censoring variable (C) is generated from a \eqn{Uniform[0,cens.rate]} distribution.
#' @param copulafam a character string specifying the family of an Archimedean copula.
#' Currently supported families are "frank", "clayton",  "amh", "gumbel", and "joe".
#' @param tau  Kendall's tau of an Archimedean copula.
#' @param params a list of marginsDist, rate1, and rate2. marginsDist is object of class "character" specifying the marginal distributions of nonterminal (T1) and terminal (T2) event times in sequence.
#' rate1 and rate2 are parameter values of their marginal distributions.
#'
#' @return A data.frame containing semi-competing risks outcomes from \code{n} subjects.
#' It is of dimension \eqn{n\times 5} with the columns named as T1, T2, event1, event2 and tr. \cr
#' \describe{
#'   \item{T1}{a vector of the observed nonterminal event times for the n individuals, i.e.,
#'    \eqn{\min(\tilde{T}_1,\tilde{T}_2,C)}, where \eqn{\tilde{T}_1}, and \eqn{\tilde{T}_2}
#'    are the true nonterminal and terminal event times, respectively.}
#'   \item{T2}{a vector of the observed terminal event times for the n individuals, i.e.,
#'    \eqn{\tilde{T}_2,C)}.}
#'   \item{event1}{a vector of censoring indicators for the non-terminal event time (1=event occurred, 0=censored).}
#'   \item{event2}{a vector of censoring indicators for the terminal event time (1=event occurred, 0=censored).}
#'   \item{tr}{treatment indicator.}
#' }
#' @importFrom copula archmCopula tau iTau frankCopula gumbelCopula claytonCopula joeCopula amhCopula
#' @importFrom copula copAMH copClayton copFrank rstable1 copJoe rCopula rMvdc mvdc
#' @importFrom acopula nderive
#'
#' @export
#' @seealso \code{\link{scrassonp}}, \code{\link{scrsurv}}



simSCRtr<- function (n  = 100,  K = 2, cens.rate = 5,copulafam="frank",tau=0.5,
                     params=list(marginsDist = rep("exp",2),rate1=1,rate2=1))
{

  alpha<- sapply(tau, Calitau,copulafam = copulafam)


  tr<- rbinom(n, size = K-1, prob = 0.5)
  biSurvTimes<- matrix(NA,nrow = n,ncol = 2)

  utr<- 0:(K-1)
  if(length(alpha)==1){alpha<- rep(alpha,K)}
  if(length(params$rate1)==1){params$rate1<- rep(params$rate1,K)}
  if(length(params$rate2)==1){params$rate2<- rep(params$rate2,K)}
  if(is.vector(params$marginsDist)){
    if(length(params$marginsDist)==2){
      params$marginsDist<- rep(list(params$marginsDist),3)
    }else{
      stop("The length of marginsDist in params shall be 2.")
    }
  }else{
    if(is.list(params$marginsDist)){
      if(length(params$marginsDist)!=K){
        stop("The length of list marginsDist in params shall equal K.")
      }else{
        if(any(sapply(params$marginsDist, length)!=2)){
          stop("The length of each element of marginsDist in params shall be 2.")
        }
    }
    }else{
      stop("The form of marginsDist in params is not correct.")
    }
  }
  for (u in utr) {
    nsim<- sum(tr==u)
    au<- alpha[utr==u]
    rate1<- params$rate1[utr==u]
    rate2<- params$rate2[utr==u]
    mdist<- params$marginsDist[utr==u][[1]]
    tsim<- copula::rMvdc(nsim,copula::mvdc(copula =copula::archmCopula(family= copulafam,
                                                                    param = au,
                                                                    dim = 2),
                                        margins = mdist,
                                        paramMargins = c(list(list(rate = rate1)),
                                                         list(list(rate = rate2)))))
    biSurvTimes[tr==u,]<- tsim
  }

  T1<- biSurvTimes[,1]
  T2<- biSurvTimes[,2]

  C <- stats::runif(n,min = 0 ,max = cens.rate)

  T2obs<- pmin(T2,C)
  T1obs<- pmin(T1,T2obs)
  delta1<- 1*(T1<= T2obs)
  delta2<- 1*(T2<= C)


  data<- data.frame(T1= T1obs,T2= T2obs,event1= delta1, event2= delta2,tr=tr)
  return(data)

}
