#################################################################
##   CopulaSCR R package by Tonghui Yu Copyright (C) 2024
##
## The below function is used for association and marginal analysis
## for semi-competing risks data with multiple intermediate event times
#################################################################

#' Association Analysis for Semi-Competing Risks Data with Multiple Intermediate Event Times
#'
#' Fits a copula-based model for semi-competing risks data with multiple
#' intermediate event times. The function first estimates the marginal survival
#' curves and pairwise association parameters for each intermediate event and the
#' terminal event, and then estimates the dependence parameter among multiple
#' intermediate event times.
#' @usage mscr(mT1,T2,mevent1,event2, 
#'             copulafam=c("frank","clayton","joe","gumbel","amh"), 
#'             nsim=500,ncore= 1,tol=0.01,a=NULL,b = NULL, 
#'             positive=TRUE, msurv.method= "JFKC",lower=0.02,upper=0.96)
#' @param mT1 A matrix or data frame of observed intermediate event times.
#'   Each column corresponds to one intermediate event type, and the column names,
#'   if present, are treated as event names.
#' @param T2 A numeric vector of observed terminal event times. Its length must
#'   equal the number of rows of \code{mT1}.
#' @param mevent1 A matrix or data frame of censoring indicators for intermediate
#'   event times. Its dimensions must match those of \code{mT1}.
#' @param event2 A numeric or logical vector of censoring indicators for the
#'   terminal event time, where 1 indicates the terminal event is observed and
#'   0 indicates censoring.
#' @param copulafam A character string specifying the Archimedean copula family.
#'   Supported families include \code{"frank"}, \code{"clayton"},
#'   \code{"joe"}, \code{"gumbel"}, and \code{"amh"}.
#' @param nsim Monte Carlo sample size used in the numerical approximation of the
#'   pseudo-likelihood.
#' @param ncore Number of CPU cores used for parallel computation.
#' @param tol Desired numerical tolerance in the estimation procedure.
#' @param a A positive constant used to dampen the weight function in the
#'   marginal association estimation. If \code{NULL}, it is chosen internally.
#' @param b A positive constant used to dampen the weight function in the
#'   marginal association estimation. If \code{NULL}, it is chosen internally.
#' @param positive Logical; whether to impose a positive constraint on the
#'   marginal Kendall's tau parameter in the pairwise association analysis.
#' @param msurv.method A character string specifying the method used to estimate
#'   the marginal survival curves for intermediate event times. Supported methods
#'   include \code{"JFKC"}, \code{"LRA"}, \code{"FJC"}, and \code{"Ker"}.
#' @param lower Lower bound of Kendall's tau corresponding to the copula
#'   parameter \eqn{\alpha} for the joint distribution of multiple intermediate
#'   event times.
#' @param upper Upper bound of Kendall's tau corresponding to the copula
#'   parameter \eqn{\alpha}.
#'
#' @return An object of class \code{"mscr"} representing the fitted model.
#'   The returned object is a list containing at least the following components:
#'   \describe{
#'     \item{alpha}{Estimated copula parameter for the dependence among multiple
#'       intermediate event times.}
#'     \item{logLik.alpha}{Pseudo-log-likelihood evaluated at \code{alpha}.}
#'     \item{alpha.int}{Lower and upper bounds used in the search interval for
#'       \code{alpha}.}
#'     \item{call}{The matched function call.}
#'     \item{copulafam}{The specified Archimedean copula family.}
#'     \item{theta}{Estimated copula parameters for the association between each
#'       intermediate event and the terminal event.}
#'     \item{tau.alpha}{Estimated Kendall's tau corresponding to \code{alpha}.}
#'     \item{tau.theta}{Estimated Kendall's tau values corresponding to
#'       \code{theta}.}
#'     \item{S_D}{Estimated survival function of the terminal event time.}
#'     \item{S_k}{A list of estimated survival functions for the intermediate
#'       event times.}
#'     \item{mar.fits}{A list of fitted \code{scrsurv} objects for the marginal
#'       analyses.}
#'   }
#'
#' @seealso \code{\link{plot.mscr}}, \code{\link{marginal_fit}}, \code{\link{terminal_survival}}
#'
#' @export
#'
#' @references Yu, Tonghui, and Xiang, Liming (2026).
#' Learning association from multiple intermediate events for dynamic prediction of survival:
#' an application to cardiovascular disease prognosis.
#' \emph{Biometrics}, 82(2), ujag087.
#' \doi{10.1093/biomtc/ujag087}
#'
#' @references Yu, Tonghui, Zhang, Binhui, Xiang, Liming, Ma, Jianwei, and Chen, Chixiang. CopulaSCR: An R Package for the Analysis of Semi-competing Risks Data Using Copula-based Models.
#' \emph{Working paper.}
#'
#' @examples
#' \donttest{
#' set.seed(12345)
#' data <- simSCRmul(
#'   n = 100, copulafam = "frank", K = 3,
#'   tau.alpha = 0.5, tau.theta = 0.5
#' )
#'
#' fit <- mscr(
#'   mT1 = data$T1,
#'   T2 = data$T2,
#'   mevent1 = data$event1,
#'   event2 = data$event2,
#'   copulafam = "frank",
#'   ncore = 1
#' )
#'
#' print(fit)
#' plot(marginal_fit(fit, event = 1))
#' plot(fit, type = "msurv")
#' plot(fit, type = "data")
#' }
#'
#'
mscr<- function(mT1,T2,mevent1,event2,
                copulafam=c("frank","clayton","joe","gumbel","amh"),
                nsim=500,ncore= 1,tol=0.01,a=NULL,b = NULL,
                positive=TRUE,msurv.method= "JFKC",lower=0.02,upper=0.96){
  Call<- match.call()
  copulafam<- match.arg(copulafam)
  if (!copulafam %in% c("frank", "clayton", "gumbel")) {
    warning("Copula families other than 'frank', 'clayton', and 'gumbel' may be numerically unstable in this implementation.",
            call. = FALSE)
  }
  indx<- match(c("mT1","T2","mevent1","event2"), names(Call), nomatch=0)

  if (indx[1]==0) stop("a mT1 argument is required")
  if (indx[2]==0) stop("a T2 argument is required")
  if (indx[3]==0) stop("a mevent1 argument is required")
  if (indx[4]==0) stop("a event2 argument is required")

  if (!(is.matrix(mT1) | is.data.frame(mT1))) {
    stop("mT1 object must be a matrix or a data.frame")
  }

  if (!(is.matrix(mevent1) | is.data.frame(mevent1))) {
    stop("mevent1 object must be a matrix or a data.frame.")
  }

  if (!(is.matrix(T2) | is.vector(T2))) {
    stop("T2 object must be a matrix or a vector.")
  }

  if (!(is.matrix(event2) | is.vector(event2))) {
    stop("event2 object must be a matrix or a vector.")
  }

  T2<- as.vector(T2)
  event2<- as.vector(event2)

  if(ncol(mT1)!=ncol(mevent1)){
    stop("mT1 and mevent1 shall have same number of columns and rows.")
  }

  if(nrow(mT1)!=nrow(mevent1)){
    stop("mT1 and mevent1 shall have same number of columns and rows.")
  }

  if(nrow(mT1)!=length(T2)){
    stop("The length of T2 shall equal the number of rows of mT1.")
  }

  if(length(event2)!=length(T2)){
    stop("T2 and event2 shall have same length.")
  }

  K<- ncol(mT1)

  if(K==1){
    stop("The number of columns in mT1 and mevent1 shall >1.")
  }

  t1fits<- list()

  if(ncore>1&!(K<=3&sum(1-event2)<10)){
    Cores<- parallel::makeCluster(min(parallel::detectCores(),ncore,2^K))
    doParallel::registerDoParallel(Cores)
    parallel<- TRUE
  }else{
    parallel<- FALSE
  }


  if(isTRUE(parallel)){
    t1fits<- foreach(k = 1:ncol(mT1),
                     .export = c("scrassonp","scrassonp0","assonpfit","Caltau","scrnpGEE",
                                 "npCHRfit","prepformula","prepData0","prepData",
                                 "copulabd_pos","copulabd"))%dopar%{
                                   sdata<- data.frame(T1 = mT1[,k],T2=T2,event1=mevent1[,k],event2=event2)
                                   fitasso<- scrassonp(t1.formula=Surv(T1, event1) ~ 1,
                                                       t2.formula=Surv(T2, event2) ~ 1,
                                                       data = sdata,copulafam=copulafam, B=0,
                                                       se = FALSE,model=TRUE,a = a,b = b,tol=1e-6,positive = positive)
                                   scrsurv(fit = fitasso, method= msurv.method,
                                           surv2km=TRUE,conf.int =FALSE)
                                 }
  }else{
    for(k in 1:ncol(mT1)){
      sdata<- data.frame(T1 = mT1[,k],T2=T2,event1=mevent1[,k],event2=event2)
      fitasso<- scrassonp(t1.formula=Surv(T1, event1) ~ 1,
                          t2.formula=Surv(T2, event2) ~ 1,
                          data = sdata,copulafam=copulafam,
                          se = FALSE,model=TRUE,a = a,b = b,tol=1e-6,positive = positive)

      # cat("tau.est =", fitasso$tau,"\n")
      t1fits[[k]]<- scrsurv(fit = fitasso, method= msurv.method,
                            surv2km=TRUE,conf.int =FALSE)
      rm(fitasso,sdata)
    }
  }


  # S_D<- approxfun(x= c(0,t1fits[[1]]$t2.surv$time),y=c(1,t1fits[[1]]$t2.surv$surv),
  #                 yleft = 1,yright = min(t1fits[[1]]$t2.surv$surv),method = "linear",
  #                 f = 0,ties = mean,na.rm = TRUE)
  # S_k<- lapply(1:K,function(k)
  #   approxfun(x=c(0,t1fits[[k]]$t1.surv$time),y=c(1,t1fits[[k]]$t1.surv$surv),
  #             yleft = 1,yright = min(t1fits[[k]]$t1.surv$surv),method = "linear",
  #             f = 0,ties = mean,na.rm = TRUE))


  S_D<- stepfun(x= t1fits[[1]]$t2.surv$time,y=c(1,t1fits[[1]]$t2.surv$surv),
                right = FALSE,ties = max)
  S_k<- lapply(1:K,function(k)
    stepfun(x=t1fits[[k]]$t1.surv$time,y=c(1,t1fits[[k]]$t1.surv$surv),
            right = FALSE,ties = max))

  names(t1fits)<- names(S_k)<- colnames(mT1)

  theta<- sapply(t1fits,function(x)x$copulaparam)
  Sdev_k<- lapply(1:K,function(k)
    survnumdiff(t1fits[[k]]$t1.surv$time,t1fits[[k]]$t1.surv$surv))
  Sdev_D<- survnumdiff(t1fits[[1]]$t2.surv$time,t1fits[[1]]$t2.surv$surv)
  # Sdev_k<- lapply(1:K,function(k)
  #   stepfun(x=t1fits[[k]]$t1.surv$time,
  #           y=c(1,diff(c(1,t1fits[[k]]$t1.surv$surv))),
  #           right = FALSE))
  # Sdev_D<- stepfun(x=t1fits[[1]]$t2.surv$time,
  #                  y=c(1,diff(c(1,t1fits[[1]]$t2.surv$surv))),
  #                  right = FALSE)

  alphafit<- mscrassofit(mT1 = mT1,T2 = T2,mevent1 = mevent1,event2 = event2,
                         S_D = S_D,S_k = S_k,Sdev_D = Sdev_D,Sdev_k = Sdev_k,
                         theta = theta,copulafam = copulafam,nsim=nsim,
                         tol=1e-6,lower=lower,upper=upper)

  if(isTRUE(parallel)){
    parallel::stopCluster(Cores)
  }

  result<- alphafit
  result$alpha.int<- c(lower,upper)
  result$call<- Call
  result$copulafam<- copulafam
  result$theta<- theta
  result$tau.alpha<- Caltau(copulafam,alphafit$alpha)
  result$tau.theta<- sapply(theta,Caltau,copulafam =copulafam)
  result$S_D<- S_D
  result$S_k<- S_k
  result$mar.fits<- t1fits
  result$t1obs<- as.matrix(mT1)
  result$t1obs[as.matrix(mevent1) == 0]<- NA_real_
  class(result)<- "mscr"
  return(result)

}

#' @noRd
mscrassofit<- function(mT1, T2, mevent1, event2,
                       S_D, S_k, Sdev_D, Sdev_k,
                       theta, copulafam,
                       nsim = 100, tol = 1e-6,
                       lower = 0.02, upper = 0.98) {

  n<- length(T2)
  K<- ncol(mT1)
  di<- apply(mevent1, 1, sum, na.rm = TRUE)

  SD.Y<- sapply(T2, S_D)
  Sdev.Y<- sapply(T2, Sdev_D)

  Sk.T<- lapply(1:K, function(k) sapply(mT1[, k], S_k[[k]]))
  Sdevk.t<- lapply(1:K, function(k) sapply(mT1[, k], Sdev_k[[k]]))

  S_k_D<- lapply(1:K, function(k) {
    sapply(1:n, function(i) {
      Calcondsurv_D(
        surv1 = Sk.T[[k]][i],
        surv2 = SD.Y[i],
        copulafam = copulafam,
        copulaparam = theta[k]
      )
    })
  })

  Sdev_k_D <- lapply(1:K, function(k)
    sapply(1:n, function(i)
      dev_Copula(
        copulafam = copulafam,
        param = theta[k],
        p1 = Sk.T[[k]][i],
        p2 = SD.Y[i],
        mode = "12"
      )
    ) * Sdevk.t[[k]] * (mT1[, k] <= T2)
  )

  jps2<- sort(unique(c(0, T2, 2 * T2)), decreasing = FALSE)
  jps1<- lapply(1:K, function(k) sort(unique(c(0, mT1[, k]))))

  SD.t<- sapply(jps2, S_D)
  keep<- c(1, which(diff(SD.t) < 0) + 1)
  jps2<- jps2[keep]
  SD.t<- SD.t[keep]

  pos1<- which(event2 == 1)
  pos0<- which(event2 == 0)

  interval<- copulabd_mT1(
    copulafam = copulafam,
    lower = upper,
    upper = lower
  )

  prep_L1<- precompute_Calplik1(
    K = K,
    mevent1 = mevent1,
    di = di,
    Sdev.Y = Sdev.Y,
    S_k_D = S_k_D,
    Sdev_k_D = Sdev_k_D
  )

  if (length(pos0) > 0L) {
    prep_L2<- precompute_Calplik2(
      mT1 = mT1,
      mevent1 = mevent1,
      T2 = T2,
      pos0 = pos0,
      jps2 = jps2,
      SD.t = SD.t,
      S_k = S_k,
      Sdev_k = Sdev_k,
      copulafam = copulafam,
      theta = theta
    )
  } else {
    prep_L2<- list()
  }

  optfit<- optimize(
    f = function(alpha) {
      -Calplik(
        alpha = alpha,
        n = n,
        pos0 = pos0,
        pos1 = pos1,
        copulafam = copulafam,
        prep_L1 = prep_L1,
        prep_L2 = prep_L2,
        nsim = nsim
      )
    },
    interval = interval,
    tol = tol
  )

  alpha<- optfit$minimum
  logLik<- -optfit$objective

  list(alpha = alpha, logLik.alpha = logLik)
}

#' @noRd
copulabd_mT1<- function(copulafam,lower=0.02,upper=0.98){
  # lb<- switch(copulafam,
  #             amh = -1,
  #             clayton = -1,
  #             frank = 0.001,
  #             gumbel = 1,
  #             joe = 1)
  # ub<- switch(copulafam,
  #             amh = 1,
  #             clayton = 400,
  #             frank = 400,
  #             gumbel = 400,
  #             joe = 400)
  # return(c(lb,ub))
  return(c(Calitau(copulafam,lower),Calitau(copulafam,upper)))
}

#' @noRd
Calplik<- function(alpha, n, pos0, pos1, 
                   copulafam, prep_L1, prep_L2,
                   nsim = 500) {

  xsim<- lpxsim(n = nsim, copulafam = copulafam, param = alpha)

  L<- rep(1, n)

  if (length(pos1) > 0L) {
    L1<- Calplik1(
      alpha = alpha,
      prep = prep_L1,
      copulafam = copulafam,
      xsim = xsim
    )
    L[pos1]<- L1[pos1]
  }

  if (length(pos0) > 0L) {
    L[pos0]<- Calplik2(
      alpha = alpha,
      prep = prep_L2,
      copulafam = copulafam,
      xsim = xsim
    )
  }

  if (copulafam == "clayton") {
    L <- L[L > .Machine$double.xmin & L < 1.2]
  } else {
    L <- L[L < 1.2]
  }

  loglik<- log(L)
  loglik[is.infinite(loglik)]<- NA_real_

  if (sum(!is.na(loglik)) > 1L) {
    mean(loglik[!is.na(loglik)])
  } else {
    -.Machine$double.xmax
  }
}

#' @noRd
Calplik1<- function(alpha, prep, copulafam, xsim) {
  n<- prep$n
  K<- prep$K

  phi_mat<- matrix(0, nrow = n, ncol = K)
  dev_mat<- matrix(1, nrow = n, ncol = K)

  for (k in seq_len(K)) {
    p_k<- prep$S_k_D_mat[, k]
    phi_mat[, k]<- archmCopulaLink(
      copulafam = copulafam,
      param = alpha,
      p = p_k
    )

    dev_mat[, k]<- archmCopulaLink_dev(
      copulafam = copulafam,
      param = alpha,
      js = p_k
    ) * prep$Sdev_k_D_mat[, k]
  }

  phiallK<- rowSums(phi_mat, na.rm = TRUE)

  dev_mat[prep$mevent0_mask]<- 1
  phiprod<- apply(dev_mat, 1, prod, na.rm = TRUE)

  psi_d<- numeric(n)

  for (d in prep$unique_di) {
    idx<- prep$di_groups[[as.character(d)]]
    if (length(idx) == 0L) next

    weight_d<- (-xsim)^d
    mat_d<- exp(-outer(xsim, phiallK[idx]))
    psi_d[idx]<- colMeans(mat_d * matrix(weight_d, nrow = length(xsim), ncol = length(idx)),
                          na.rm = TRUE)
  }

  L1<- psi_d * phiprod * prep$Sdev.Y * prep$sign_term
  L1[is.na(L1)]<- 0
  L1
}

#' @noRd
precompute_Calplik1<- function(K, mevent1, di, Sdev.Y, S_k_D, Sdev_k_D) {
  n<- length(di)

  S_k_D_mat<- do.call(cbind, S_k_D)
  Sdev_k_D_mat<- do.call(cbind, Sdev_k_D)

  if (!is.matrix(S_k_D_mat)) {
    S_k_D_mat<- matrix(S_k_D_mat, nrow = n, ncol = K)
  }
  if (!is.matrix(Sdev_k_D_mat)) {
    Sdev_k_D_mat<- matrix(Sdev_k_D_mat, nrow = n, ncol = K)
  }

  unique_di<- sort(unique(di))
  di_groups<- lapply(unique_di, function(d) which(di == d))
  names(di_groups)<- as.character(unique_di)

  list(
    n = n,
    K = K,
    di = di,
    unique_di = unique_di,
    di_groups = di_groups,
    mevent0_mask = (mevent1 == 0),
    Sdev.Y = Sdev.Y,
    sign_term = vapply(1 + di, powm1, numeric(1)),
    S_k_D_mat = S_k_D_mat,
    Sdev_k_D_mat = Sdev_k_D_mat
  )
}

#' @noRd
Calplik2<- function(alpha, prep, copulafam, xsim) {
  L2<- numeric(length(prep))
  for (kk in seq_along(prep)) {
    prekk<- prep[[kk]]
    if (is.null(prekk$ny) || prekk$ny == 0L) {
      L2[kk]<- 0
      next
    }
    loc1_alpha<- loc1_alpha(alpha = alpha, prekk = prekk, copulafam = copulafam)
    vals<- vapply(prekk$loc0_sets, function(loc0vec) {
      L2_Jint(alpha = alpha, loc0 = loc0vec,
              prekk = prekk, loc1_alpha = loc1_alpha,
              copulafam = copulafam, xsim = xsim)
    }, numeric(1))
    L2[kk]<- sum(vals, na.rm = TRUE)
  }
  L2
}

#' @noRd
precompute_Calplik2<- function(mT1, mevent1, T2, pos0,
                               jps2, SD.t, S_k, Sdev_k,
                               copulafam, theta) {
  prep<- vector("list", length(pos0))
  for (kk in seq_along(pos0)) {
    common<- prep_subject_common(
      kk = kk, pos0 = pos0, mT1 = mT1, mevent1 = mevent1, T2 = T2,
      jps2 = jps2, SD.t = SD.t, S_k = S_k, Sdev_k = Sdev_k
    )
    prep[[kk]]<- precompute_common(common, theta = theta, copulafam = copulafam)
  }
  prep
}

#' @noRd
loc1_alpha<- function(alpha, prekk, copulafam) {
  dd<- prekk$dd
  if (dd == 0L) return(list(dd = 0L, phi1 = NULL, phiprod = NULL))

  phi1<- archmCopulaLink(copulafam = copulafam, param = alpha, p = prekk$D1_mode2[[1]])
  if (dd > 1L) for (j in 2:dd) {
    phi1<- phi1 + archmCopulaLink(copulafam = copulafam, param = alpha, p = prekk$D1_mode2[[j]])
  }

  phiprod<- archmCopulaLink_dev(copulafam = copulafam, param = alpha, js = prekk$D1_mode2[[1]]) * prekk$D1_mode12[[1]]
  if (dd > 1L) for (j in 2:dd) {
    phiprod<- phiprod *
      (archmCopulaLink_dev(copulafam = copulafam, param = alpha, js = prekk$D1_mode2[[j]]) *
         prekk$D1_mode12[[j]])
  }

  list(dd = dd, phi1 = phi1, phiprod = phiprod)
}

#' @noRd
L2_Jint<- function(alpha, loc0, prekk, loc1_alpha, copulafam, xsim) {
  ny<- prekk$ny
  if (is.null(ny) || ny == 0L) return(0)

  SD.t2<- prekk$SD.t2
  K<- prekk$K

  loc0<- as.integer(loc0)
  loc1<- prekk$loc1
  dd<- loc1_alpha$dd

  loc00<- setdiff(seq_len(K), union(loc0, loc1))

  if (length(loc00) > 0L) {
    phi00<- archmCopulaLink(copulafam = copulafam, param = alpha, p = prekk$D_tt[loc00[1], ])
    if (length(loc00) > 1L) for (ii in 2:length(loc00)) {
      d<- loc00[ii]
      phi00<- phi00 + archmCopulaLink(copulafam = copulafam, param = alpha, p = prekk$D_tt[d, ])
    }
  } else {
    phi00<- rep(0, ny)
  }

  if (dd == 0L && length(loc0) == 0L) {
    psi_d<- archmCopulaLink_inv(copulafam = copulafam, param = alpha, y = phi00)
    return(sum(diff(c(SD.t2, 0)) * psi_d, na.rm = TRUE) * (-1))
  }

  nx<- length(xsim)

  if (length(loc0) == 0L && dd != 0L) {
    phi1<- loc1_alpha$phi1
    phiprod<- loc1_alpha$phiprod

    baseExp<- exp(-outer(xsim, phi1 + phi00))
    weight <- (-xsim)^dd
    mat    <- baseExp * matrix(weight, nrow = nx, ncol = ny)

    psi_d<- colMeans(mat, na.rm = TRUE)
    Js<- sum(diff(c(SD.t2, 0)) * psi_d * phiprod, na.rm = TRUE) * powm1(dd + 1)

    if (!is.na(Js) && Js > 1) {
      med<- apply(mat, 2, median, na.rm = TRUE)
      Js<- sum(diff(c(SD.t2, 0)) * med * phiprod, na.rm = TRUE) * powm1(dd + 1)
    }
    return(Js)
  }

  prodMat<- matrix(1.0, nrow = nx, ncol = ny)
  for (d in loc0) {
    a<- archmCopulaLink(copulafam = copulafam, param = alpha, p = prekk$D_tt[d, ])
    b<- archmCopulaLink(copulafam = copulafam, param = alpha, p = prekk$D_Yt[d, ])
    prodMat<- prodMat * (exp(-outer(xsim, a)) - exp(-outer(xsim, b)))
  }

  if (dd == 0L) {
    baseExp<- exp(-outer(xsim, phi00))
    phi.t<- colMeans(baseExp * prodMat, na.rm = TRUE)
    return(sum(diff(c(SD.t2, 0)) * phi.t, na.rm = TRUE) * powm1(1 + length(loc0)))
  } else {
    phi1<- loc1_alpha$phi1
    phiprod<- loc1_alpha$phiprod

    baseExp<- exp(-outer(xsim, phi1 + phi00))
    weight <- (-xsim)^dd
    matAll <- baseExp * prodMat * matrix(weight, nrow = nx, ncol = ny)
    phi.t  <- colMeans(matAll, na.rm = TRUE)

    return(sum(diff(c(SD.t2, 0)) * phi.t * phiprod, na.rm = TRUE) *
             powm1(dd + 1 + length(loc0)))
  }
}

#' @noRd
prep_subject_common <- function(kk, pos0, mT1, mevent1, T2, jps2, SD.t, S_k, Sdev_k) {
  t1obs <- mT1[pos0[kk], ]
  t2obs <- T2[pos0[kk]]
  K <- length(t1obs)

  loc0_full <- which(mevent1[pos0[kk], ] == 0)
  loc1      <- which(mevent1[pos0[kk], ] == 1)

  idx2 <- which(jps2 >= t2obs)
  ny <- length(idx2)
  if (ny == 0L) return(list(ny = 0L))

  jps2_2 <- jps2[idx2]
  SD.t2  <- SD.t[idx2]

  Sk_grid <- matrix(NA_real_, nrow = K, ncol = ny)
  for (d in 1:K) Sk_grid[d, ] <- S_k[[d]](jps2_2)

  Sk_t2   <- numeric(K)
  Sk_t1   <- numeric(K)
  Sdev_t1 <- numeric(K)
  for (d in 1:K) {
    Sk_t2[d]   <- S_k[[d]](t2obs)
    Sk_t1[d]   <- S_k[[d]](t1obs[d])
    Sdev_t1[d] <- Sdev_k[[d]](t1obs[d])
  }

  # subsets of loc0 (include empty)
  loc0_sets <- list(integer(0))
  if (length(loc0_full) > 1) loc0_sets <- c(loc0_sets, get_subsets(loc0_full))
  else if (length(loc0_full) == 1) loc0_sets <- c(loc0_sets, list(loc0_full))

  list(
    ny = ny, K = K,
    t1obs = t1obs, t2obs = t2obs,
    loc1 = as.integer(loc1), dd = length(loc1),
    loc0_sets = loc0_sets,
    SD.t2 = SD.t2,
    Sk_grid = Sk_grid,
    Sk_t2 = Sk_t2,
    Sk_t1 = Sk_t1,
    Sdev_t1 = Sdev_t1
  )
}

#' @noRd
precompute_common<- function(common, theta, copulafam) {
  ny<- common$ny
  if (is.null(ny) || ny == 0L) return(list(ny = 0L))

  K<- common$K
  SD.t2<- common$SD.t2
  Sk_grid<- common$Sk_grid
  Sk_t2<- common$Sk_t2
  Sk_t1<- common$Sk_t1
  Sdev_t1<- common$Sdev_t1
  loc1<- common$loc1
  dd<- common$dd

  # alpha-free derivatives
  D_tt<- matrix(NA_real_, nrow = K, ncol = ny)
  D_Yt<- matrix(NA_real_, nrow = K, ncol = ny)
  for (d in 1:K) {
    D_tt[d, ]<- dev_Copula(copulafam = copulafam, param = theta[d],
                           p1 = Sk_grid[d, ], p2 = SD.t2, mode = "2")
    D_Yt[d, ]<- dev_Copula(copulafam = copulafam, param = theta[d],
                           p1 = Sk_t2[d], p2 = SD.t2, mode = "2")
  }

  if (dd > 0L) {
    D1_mode2 <- vector("list", dd)
    D1_mode12<- vector("list", dd)
    for (j in 1:dd) {
      d<- loc1[j]
      D1_mode2[[j]]<- dev_Copula(copulafam = copulafam, param = theta[d],
                                 p1 = Sk_t1[d], p2 = SD.t2, mode = "2")
      D1_mode12[[j]]<- dev_Copula(copulafam = copulafam, param = theta[d],
                                  p1 = Sk_t1[d], p2 = SD.t2, mode = "12") * Sdev_t1[d]
    }
  } else {
    D1_mode2<- list()
    D1_mode12<- list()
  }

  # merge: keep common + derivative caches
  within(common, {
    D_tt<- D_tt
    D_Yt<- D_Yt
    D1_mode2<- D1_mode2
    D1_mode12<- D1_mode12
  })
}

#' @noRd
lpxsim<- function(n,copulafam,param){
  if(copulafam=="amh"){
    xsim<- copula::copAMH@V0(n = n, theta = param) # geometric
  }
  if(copulafam=="clayton"){
    xsim<-copula::copClayton@V0(n = n, theta = param) # gamma
  }
  if(copulafam=="frank"){
    xsim<-copula::copFrank@V0(n = n, theta = param) #logarithmic series
  }
  if(copulafam=="gumbel"){
    xsim<-copula::rstable1(n = n,alpha = 1/param,beta = 1,
                           gamma = (cos(pi/2/param))^{param},delta = 0,pm = 1)
    # positive stable
  }
  if(copulafam=="joe"){
    xsim<-copula::copJoe@V0(n = n, theta = param) # Sibuya latent (power series)
  }

  return(xsim)
}



# psialpha_dev<- function(d,copulafam,param){
#   # psialpha_dev<- function(u,d,copulafam,param){
#   # if(copulafam=="amh"){
#   #   phi_inv<-  "(param-1)/(param - exp(y))"
#   # }
#   # if(copulafam=="clayton"){
#   #   phi_inv<- "(param*y+1)^{-1/param}"
#   # }
#   # if(copulafam=="frank"){
#   #   phi_inv<- "-1/param*log(1+(exp(-param)-1)*exp(y))"
#   # }
#   # if(copulafam=="gumbel"){
#   #   phi_inv<- "exp(-y^(1/param))"
#   # }
#   # if(copulafam=="joe"){
#   #   phi_inv<-"1- (1-exp(-y))^(1/param)"
#   # }
#   # calculus::derivative(f = phi_inv, var = c(y = u),
#   #                      order = d,params = list(param=param))
#
#
#   if(copulafam=="amh"){
#     phi_inv<- expression((param-1)/(param - exp(y)))
#   }
#   if(copulafam=="clayton"){
#     phi_inv<- expression((param*y+1)^{-1/param})
#   }
#   if(copulafam=="frank"){
#     phi_inv<- expression(-1/param*log(1+(exp(-param)-1)*exp(y)))
#   }
#   if(copulafam=="gumbel"){
#     phi_inv<- expression(exp(-y^(1/param)))
#   }
#   if(copulafam=="joe"){
#     phi_inv<-expression(1- (1-exp(-y))^(1/param))
#   }
#
#   df<-DD(phi_inv,"y",d)
#
#
#   df
# }
# DD <- function(expr, name, order = 1,...) {
#   if(order < 1) stop("'order' must be >= 1")
#   if(order == 1) D(expr, name,...)
#   else DD(D(expr, name,...), name, order - 1,...)
# }
#' @noRd
get_subsets<- function(vec){
  vec<- as.numeric(unique(vec))
  unlist(lapply(1:length(vec),
                combn,
                x = vec,
                simplify = FALSE),
         recursive = FALSE)

}
#' @noRd
Calcondsurv_D<- function(surv1,surv2,copulafam,copulaparam){
  if(length(surv1)>1&length(surv2)>1){
    out<- lapply(surv2, function(xx)
      dev_Copula(copulafam = copulafam,param = copulaparam,p1 =surv1,p2 =xx,mode = "2"))
  }else{
    out<- dev_Copula(copulafam = copulafam,param = copulaparam,p1 =surv1,p2 =surv2,mode = "2")
  }
  return(out)

}

#' @noRd
survnumdiff<- function(x,y,smooth=FALSE,jumpsize=FALSE,span = 0.75,exact=TRUE){

  if(!isTRUE(exact)){

    # da<- data.frame(x=x.new,y=y.new)
    if(isTRUE(jumpsize)){
      y<- c(1,y)
      x<- c(0,x)
      # x.new<- x[!is.na(y)]
      # y.new<- y[!is.na(y)]
      ypos<- unique(y[!is.na(y)])
      ypos<- sapply(ypos,function(yy)which(y==yy)[1])
      x.new<- x[ypos]
      y.new<- y[ypos]
      sdev<-  diff(c(1,y.new),lag = 1)
    }else{
      y<- c(1,y)
      x<- c(0,x)
      ypos<- unique(y[!is.na(y)])
      ypos<- sapply(ypos,function(yy)which(y==yy)[1])
      x.new<- x[ypos]
      y.new<- y[ypos]
      sdev<-  diff(y.new,lag = 2)/diff(x.new,lag = 2)
      sdev<- c((y.new[2]-y.new[1])/(x.new[2]-x.new[1]),sdev)
      sdev<- c(sdev,(y.new[length(y.new)]-y.new[length(y.new)-1])/
                 (x.new[length(x.new)]-x.new[length(x.new)-1]))
      x.new<- x.new[!is.na(sdev)&!is.infinite(sdev)]
      sdev<- sdev[!is.na(sdev)&!is.infinite(sdev)]

    }
    if(isTRUE(smooth)){
      sdev<- predict(stats::loess(sdev~x.new,span=span),x.new)
    }

    sdev<- approxfun(x=x.new,y=sdev,method = "linear",
                     yleft = 0,yright = 0,f=0,na.rm = TRUE,ties = mean)
  }

  if(isTRUE(exact)){
    y<- c(1,y)
    x<- c(0,x)
    x.new<- x[!is.na(y)]
    y.new<- y[!is.na(y)]

    xpos<- unique(x.new)
    xpos<- sapply(xpos,function(xx)which(x.new==xx)[1])
    x.new<- x.new[xpos]
    y.new<- y.new[xpos]

    # sdev0<-  diff(c(1,y.new),lag = 1)


    # sdev<- function(xx){
    #   xpos<- match(xx,x.new,nomatch = FALSE)
    #   yy<- numeric(length = length(xx))
    #   yy[xpos!=0]<- sdev0[xpos[xpos!=0]]
    #   yy
    # }
    # sdev<- approxfun(x=x.new,y=sdev0,method = "constant",
    #                  yleft = 0,yright = 0,f=0,na.rm = TRUE,ties = mean)


    # sdev<- stepfun(x=x.new,y=c(0,diff(c(1,y.new))), right = FALSE,ties = mean )
    sdev<- approxfun(x=x.new,y=diff(c(1,y.new)), method = "constant",
                     yleft = 0,yright = 0,f=0,na.rm = TRUE,ties = mean )
  }



  return(sdev)

}

#' @noRd
powm1<- function(p) if ((p %% 2L) == 0L) 1 else -1
