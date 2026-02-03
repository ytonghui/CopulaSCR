#' Print objects in the scrsurv library
#'
#' Pretty printing of objects created with the functionality of the `scrsurv'
#' library.
#'
#'
#' @aliases print.scrsurv
#' @param x Object of class \code{scrsurv}.
#' @param \dots Not used.
#' @seealso \code{\link{scrsurv}}, \code{\link{scrassonp}}, \code{\link{plot.scrsurv}}, \code{\link{print.scrassonp}}
#' @keywords survival
#' @export
#'
print.scrsurv <- function(x,...) {
  cat("\n")
  cat("Call: ")
  print(x$Call)
  cat("\n")

  message(switch(x$t1.surv$covariate.type,"nonparametric estimator",
                 "Stratified nonparametric estimator",
                 "Stone-Beran-type estimator",
                 "Stratified Stone-Beran-type estimator")," for the ",
          ifelse(x$covariate.type==1," "," conditional "),
          "survival functions of nonterminal and terminal event times in the semi-competing risks data. ")

  cat("\n")
  ##   discrete.predictors <- extract.name.from.special(grep("strata.",names(x$X),value=TRUE),pattern="strata\\.")
  ##   continuous.predictors <- extract.name.from.special(grep("NN.",names(x$X),value=TRUE),pattern="NN\\.")
  discrete.predictors <- x$t1.surv$discrete.predictors
  continuous.predictors <- x$t1.surv$continuous.predictors

  message(#"Predictor space:\n\n",
    switch(x$t1.surv$covariate.type,
           "No covariates",{
             if (length(discrete.predictors)==1){
               c("Discrete predictor variable: ", discrete.predictors)
             }else{
               c("Discrete predictor variables:\n",
                 sapply(discrete.predictors,function(x)paste("\n - ",x)))
             }},
           c("Continuous predictors: ",continuous.predictors),
           c("  Discrete predictor variables: ",
             paste(discrete.predictors,collapse=", "),
             "\nContinuous predictor variables: ",
             continuous.predictors)))


  nonterminal<- c("n"=as.character(length(x$t1.surv$model.response[,2])),
                  "events"=as.character(sum(x$t1.surv$model.response[,2])),
                  "maxtime" = round(max(x$t1.surv$time),2))
  terminal<- c("n"= as.character(length(x$t2.surv$model.response[,2])),
                  "events"=as.character(sum(x$t2.surv$model.response[,2])),
                  "maxtime" = round(max(x$t2.surv$time),2))

  print(rbind(nonterminal,terminal),quote = FALSE)


  if (!is.null(x$na.action)){
    cat("\n",
        length(x$na.action),
        ifelse(length(x$na.action)==1,
               " observation",
               " observations")," deleted due to missing values.\n",sep="")
  }
  invisible(x)

}


#' Print objects in the scrassonp library
#'
#' Pretty printing of objects created with the functionality of the `scrassonp'
#' library.
#'
#'
#' @aliases print.scrassonp
#' @param x Object of class \code{scrassonp}.
#' @param \dots Not used.
#' @seealso \code{\link{scrassonp}}, \code{\link{scrsurv}}
#' @keywords survival
#' @export
#'
print.scrassonp <- function(x,...) {
  cat("\n")
  cat("Call: ")
  print(x$Call)
  cat("\n")

  if(x$t1data$covariate.type==1){
    cat("Association analysis of semi-competing risks data.\n")
  }
  if(x$t1data$covariate.type>1){
    cat("Stratified association analysis of semi-competing risks data.\n")
  }

  ##   discrete.predictors <- extract.name.from.special(grep("strata.",names(x$X),value=TRUE),pattern="strata\\.")
  ##   continuous.predictors <- extract.name.from.special(grep("NN.",names(x$X),value=TRUE),pattern="NN\\.")
  discrete.predictors <- x$t1data$discrete.predictors
  continuous.predictors <- x$t1data$continuous.predictors

  message(#"Predictor space:\n\n",
    switch(x$t1data$covariate.type,
           "",{
             if (length(discrete.predictors)==1){
               c("Discrete predictor variable: ", discrete.predictors)
             }else{
               c("Discrete predictor variables:\n",
                 sapply(discrete.predictors,function(x)paste("\n - ",x)))
             }},
           c("Continuous predictors: ",continuous.predictors),
           c("  Discrete predictor variables: ",
             paste(discrete.predictors,collapse=", "),
             "\nContinuous predictor variables: ",
             continuous.predictors)))

 ##???
  cat("Estimation of Kendall's tau in the", toupper(x$copulafam),
      "copula structure:",
      round(x$tau,3),"\n")

  if(isTRUE(x$se)){
    cat("Standard error:", round(x$tau.se,3),"\n")
  }


  if (!is.null(x$t1data$na.action)){
    cat("\n",
        length(x$t1data$na.action),
        ifelse(length(x$t1data$na.action)==1,
               " observation",
               " observations")," deleted due to missing values.\n",sep="")
  }
  invisible(x)

}


#' Print objects in the mscr library
#'
#' Pretty printing of objects created with the functionality of the `scrassonp'
#' library.
#'
#'
#' @aliases print.mscr
#' @param x Object of class \code{mscr}.
#' @param \dots Not used.
#' @seealso \code{\link{mscr}}, \code{\link{scrsurv}}
#' @keywords survival
#' @export
#'
print.mscr <- function(x,...) {
  cat("\n")
  cat("Call: ")
  print(x$call)
  cat("\n")
  cat(length(x$mar.fits), "intermediate events, 1 terminal event.\n")
  cat("-------\n")
  message("Association analyses")
  cat("Estimated Kendall's tau in the", toupper(x$copulafam),
      "copula structure\n")
  cat("tau.theta:\n")
  print(round(x$tau.theta,3))

  cat("tau.alpha:\n")
  print(round(x$tau.alpha,3))
  cat("-------\n")
  message("Marginal analyses (survival curves)")
  cat("intermediate event times: A list of", length(x$mar.fits), "scrsurv objects.\n")
  cat("terminal event time:\n")
  summary(x$S_D)


  invisible(x)

}


