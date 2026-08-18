#' Plot of Marginal Survival Curves for Nonterminal and Terminal Event Times
#' @param x Object of class \code{scrsurv}.
#' @param linetype 	line types. Allowed values includes i) "strata" for changing
#'   linetypes by strata (i.e. groups); ii) a numeric vector (e.g., c(1, 2)) or
#'   a character vector c("solid", "dashed").
#' @param conf.int logical value. If TRUE, plots confidence interval.
#' @param censor 	logical value. If TRUE, censors will be drawn.
#' @param xlab axis labels
#' @param ylab axis labels
#' @param align Specifies whether graphs in the grid should be horizontally ("h")
#'    or vertically ("v") placed.
#' @param ... other arguments
## #' @import survminer
## #' @import cowplot
#' @seealso \code{\link{scrsurv}}
#' @return No return value.
#' @export
"plot.scrsurv"<- function (x,linetype = "strata",
                         conf.int = FALSE,
                         censor = FALSE,xlab = "Time",
                         ylab = "Survival probability",align="h",
                         ...)
{
  if(!is.null(x$stratas)){

    strata1<- rep(x$t1.surv[["X"]][[1]],times=x$t1.surv$size.strata)
    strata2<- rep(x$t2.surv[["X"]][[1]],times=x$t2.surv$size.strata)

    object<- x$t1.surv
    if(conf.int&!is.null(object$std.err)){

    survmat1<- data.frame(time=object$time,surv=object$surv,
                          n.event = object$n.event,
                          n.censor=object$n.censor,
                          std.err= object$std.err,
                          n.risk = object$n.risk,
                          upper=object$upper,
                          lower=object$lower,
                          strata=strata1)

    object<- x$t2.surv
    survmat2<- data.frame(time=object$time,surv=object$surv,
                          n.event = object$n.event,
                          n.censor=object$n.censor,
                          std.err= object$std.err,
                          n.risk = object$n.risk,
                          upper=object$upper,
                          lower=object$lower,
                          strata=strata2)
    }else{
      survmat1<- data.frame(time=object$time,surv=object$surv,
                            n.event = object$n.event,
                            n.censor=object$n.censor,
                            n.risk = object$n.risk,
                            strata=strata1)

      object<- x$t2.surv
      survmat2<- data.frame(time=object$time,surv=object$surv,
                            n.event = object$n.event,
                            n.censor=object$n.censor,
                            n.risk = object$n.risk,
                            strata=strata2)
      conf.int<- FALSE
    }

  }else{
    if(linetype == "strata")linetype = 1

    object<- x$t1.surv
    if(conf.int&!is.null(object$std.err)){

      survmat1<- data.frame(time=object$time,surv=object$surv,
                            n.event = object$n.event,
                            n.censor=object$n.censor,
                            std.err= object$std.err,
                            n.risk = object$n.risk,
                            upper=object$upper,
                            lower=object$lower)

      object<- x$t2.surv
      survmat2<- data.frame(time=object$time,surv=object$surv,
                            n.event = object$n.event,
                            n.censor=object$n.censor,
                            std.err= object$std.err,
                            n.risk = object$n.risk,
                            upper=object$upper,
                            lower=object$lower)
    }else{
      object<- x$t1.surv
      survmat1<- data.frame(time=object$time,surv=object$surv,
                            n.event = object$n.event,
                            n.censor=object$n.censor,
                            n.risk = object$n.risk)

      object<- x$t2.surv
      survmat2<- data.frame(time=object$time,surv=object$surv,
                            n.event = object$n.event,
                            n.censor=object$n.censor,
                            n.risk = object$n.risk)
      conf.int<- FALSE
    }



  }

  survfig<-list()
  survfig[[1]]<- survminer::ggsurvplot(survmat1,
                                       font.main = 15,font.legend=8,
                                       font.x =  10, font.y = 10,
                                       font.tickslab = 8,
                                       title="nonterminal event",
                                       linetype = linetype,conf.int = conf.int,
                                       censor = censor,xlab = xlab,ylab = ylab,
                                       ...)

  survfig[[2]]<- survminer::ggsurvplot(survmat2,
                                       font.main = 15,font.legend=8,
                                       font.x =  10, font.y = 10,
                                       font.tickslab = 8,
                                       title="terminal event",
                                       linetype = linetype,conf.int = conf.int,
                                       censor = censor,xlab = xlab,ylab = ylab,
                                       ...)
  # survfig[[2]]<- survfig[[2]]+ggplot2::ggtitle("terminal event")
  if(align=="h"){
    ncol<- 2
  }else{
    ncol<- 1
  }
  cowplot::plot_grid(plotlist=survfig,ncol=ncol)

}


#' Plot an 'mscr' Object
#' @param x Object of class \code{mscr}.
#' @param type Display the data structure if \code{type = "data"}, and
#'   display the marginal survival curve of each event time if \code{type = "msurv"}.
#' @param nt.seperate logical value. If FALSE, place all survival curves of
#'   intermediate event times in a single plot
#' @param linetype 	line types. Allowed values includes i) "strata" for changing
#'   linetypes by strata (i.e. groups); ii) a numeric vector (e.g., c(1, 2)) or
#'   a character vector c("solid", "dashed").
#' @param censor 	logical value. If TRUE, censors will be drawn.
#' @param xlab axis labels
#' @param ylab axis labels
#' @param fontsize Font size, in points, for text in the plot of the data structure.
#' @param height Height of the node in the plot of the data structure, in inches.
#' @param ... other arguments
#' @seealso \code{\link{mscr}}
#' @return No return value.
#' @export
"plot.mscr"<- function (x,type = c("msurv","data"),
                      nt.seperate = FALSE, linetype = "strata",
                      censor = FALSE, xlab = "Time",
                      ylab = "Survival probability",
                      fontsize=10,height=7,...)
{
  type<- match.arg(type)
  if(type=="data"){
    mscr_dataplot(tau.alpha=x$tau.alpha,tau.theta= x$tau.theta,
                  copulafam=x$copulafam,varnames = names(x$mar.fits),
                  fontsize=  fontsize,height=height)

  }


  if(type=="msurv"){
    mscr_msurvplot(objlist=x$mar.fits,nt.seperate=nt.seperate,
                   linetype = linetype,censor =censor,xlab = xlab,ylab = ylab,...)
  }

}

#' @noRd
mscr_msurvplot<- function(objlist,nt.seperate=TRUE,linetype,censor,xlab,ylab,...){
  survfig<-list()

  if(isTRUE(nt.seperate)){
    if(linetype == "strata")linetype = 1
    for (k in 1:length(objlist)) {
      object<- objlist[[k]]$t1.surv
      survmat1<- data.frame(time=object$time,surv=object$surv,
                            n.event = object$n.event,
                            n.censor=object$n.censor,
                            n.risk = object$n.risk)
      survfig[[k]]<-
        survminer::ggsurvplot(survmat1,
                              font.main = 15,font.legend=8,
                              font.x =  10, font.y = 10,
                              font.tickslab = 8,
                              title=names(objlist)[k],
                              linetype = linetype,conf.int = FALSE,
                              censor = censor,xlab = xlab,ylab = ylab,
                              legend = "none",
                              ...)
    }

  }else{
    survmat1<- list()
    for (k in 1:length(objlist)) {
      object<- objlist[[k]]$t1.surv
      survmat1[[k]]<-
        data.frame(time=object$time,surv=object$surv,
                   n.event = object$n.event,
                   n.censor=object$n.censor,
                   n.risk = object$n.risk,
                   strata = rep(names(objlist)[k],length(object$time)))
    }

    survmat1<- do.call(rbind,survmat1)
    survfig[[1]]<-
      survminer::ggsurvplot(survmat1,
                            font.main = 15,font.legend=8,
                            font.x =  10, font.y = 10,
                            font.tickslab = 8,
                            title="nonterminal events",
                            linetype = linetype,conf.int = FALSE,
                            censor = censor,xlab = xlab,ylab = ylab,
                            legend.title = "nonterminal events",
                            ...)
    }

  object<- objlist[[1]]$t2.surv

  survmat2<- data.frame(time=object$time,surv=object$surv,
                        n.event = object$n.event,
                        n.censor=object$n.censor,
                        n.risk = object$n.risk)




  survfig[[length(survfig)+1]]<-
    survminer::ggsurvplot(survmat2,
                          font.main = 15,font.legend=8,
                          font.x =  10, font.y = 10,
                          font.tickslab = 8,
                          title="terminal event",
                          linetype = linetype,conf.int = FALSE,
                          censor = censor,xlab = xlab,ylab = ylab,
                          legend = "none",
                          ...)

  cowplot::plot_grid(plotlist=survfig)
}



#' @noRd
mscr_dataplot<- function(tau.alpha,tau.theta,copulafam,varnames,
                         fontsize=10,height=7){
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)

  Nodes <- c("Entry",varnames,"DTH")
  Edges <- list("Entry"=list(edges=c(varnames,"DTH")))
  for (k in 1:length(varnames)) {
    Edges[[varnames[k]]]<- list(edges=c(varnames[-k],"DTH"))
    #weights= c(rep(tau.alpha,length(varnames)-1),tau.theta[k])
  }

  Edges[["DTH"]]<- list(edges=NULL)

  treeobj <- graph::graphNEL(nodes=Nodes, edgeL=Edges, edgemode="directed")

  nAttrs<-list()
  nAttrs$fontsize <- rep(fontsize, length(Nodes))
  nAttrs$fixedsize<- rep(FALSE, length(Nodes))
  nAttrs$height<- rep(height, length(Nodes))
  nAttrs$color<- c("gray70",rep("black",length(Nodes)-2),"red")


  names(nAttrs$fixedsize)<- names(nAttrs$fontsize)<-
    names(nAttrs$height)<- names(nAttrs$color)<- Nodes
  Rgraphviz::plot(treeobj,nodeAttrs=nAttrs,
              main = paste0("Under ",toupper(copulafam)," copula structure"))
  graphics::par(new = TRUE, bty = "n")
  plot(0:10,0:10,type="n",yaxt="n",xaxt="n")

  text(2,8,paste0("tau.theta= \n (",paste0(round(tau.theta,2),collapse = ","),")\n",
               "\ntau.alpha=\n",round(tau.alpha,2)),
        col="red")

}


#' Plot an 'scrassonp' Object
#'
#' @param x Object of class \code{scrassonp}.
#' @param type Plot the observed bivariate survival times if
#'   \code{type = "data"}, or display the fitted association structure if
#'   \code{type = "scrassonp"}.
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param ... Other arguments passed to the plotting function.
#'
#' @seealso \code{\link{scrassonp}}
#' @return No return value.
#' @export
plot.scrassonp <- function(x,
                           type = c("data", "scrassonp"),
                           xlab = NULL,
                           ylab = NULL,
                           ...) {
  
  type <- match.arg(type)
  
  if (type == "data") {
    
    if (is.null(xlab)) xlab <- "T1"
    if (is.null(ylab)) ylab <- "T2"
    
    nstrata <- x$t1data$nstrata
    
    draw_data <- function(ind, strata.name = NULL) {
      
      plotdata <- data.frame(
        T1 = x$t1data$time[ind],
        T2 = x$t2data$time[ind],
        event1 = x$t1data$event[ind],
        event2 = x$t2data$event[ind]
      )
      
      obj <- SurvCorr::survcorr(
        formula1 = survival::Surv(T2, event2) ~ 1,
        formula2 = survival::Surv(T1, event1) ~ 1,
        data = plotdata
      )
      
      graphics::plot(
        obj,
        "times",
        xlab = xlab,
        ylab = ylab,
        main = strata.name,
        ...
      )
    }
    
    if (nstrata == 1) {
      
      draw_data(seq_along(x$t1data$time))
      
    } else {
      
      oldpar <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(oldpar))
      
      nr <- ceiling(sqrt(nstrata))
      nc <- ceiling(nstrata / nr)
      graphics::par(mfrow = c(nr, nc))
      
      size.strata <- x$t1data$size.strata
      end.strata <- cumsum(size.strata)
      start.strata <- c(1, utils::head(end.strata, -1) + 1)
      
      strata.names <- as.character(x$stratas[[1]])
      
      for (k in seq_len(nstrata)) {
        
        ind <- start.strata[k]:end.strata[k]
        
        draw_data(
          ind = ind,
          strata.name = strata.names[k]
        )
      }
    }
  }
  
  if (type == "scrassonp") {
    
    if (is.null(xlab)) xlab <- "u"
    if (is.null(ylab)) ylab <- "v"
    
    u <- seq(0.001, 0.999, length.out = 100)
    v <- seq(0.001, 0.999, length.out = 100)
    
    params <- x$copulaparam
    taus <- x$tau
    
    draw_contour <- function(param, tau, strata.name = NULL,
                             derivative = FALSE) {
      
      uv <- expand.grid(u = u, v = v)
      
      if (!derivative) {
        
        z <- matrix(
          Copulafn(
            copulafam = x$copulafam,
            param = param,
            p1 = uv$u,
            p2 = uv$v
          ),
          nrow = length(u),
          ncol = length(v)
        )
        
        main <- paste0(
          toupper(x$copulafam),
          " copula"
        )
        
      } else {
        
        z <- matrix(
          dev_Copula(
            copulafam = x$copulafam,
            param = param,
            p1 = uv$u,
            p2 = uv$v,
            mode = "12"
          ),
          nrow = length(u),
          ncol = length(v)
        )
        
        main <- paste0(
          toupper(x$copulafam),
          " copula density"
        )
      }
      
      if (!is.null(strata.name)) {
        main <- paste0(main, " - ", strata.name)
      }
      
      main <- paste0(
        main,
        "\nKendall's tau = ",
        round(tau, 3)
      )
      
      graphics::contour(
        x = u,
        y = v,
        z = z,
        xlab = xlab,
        ylab = ylab,
        main = main,
        ...
      )
    }
    
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))
    
    if (length(params) == 1) {
      
      graphics::par(mfrow = c(1, 2))
      
      draw_contour(
        param = params,
        tau = taus,
        derivative = FALSE
      )
      
      draw_contour(
        param = params,
        tau = taus,
        derivative = TRUE
      )
      
    } else {
      
      graphics::par(mfrow = c(length(params), 2))
      
      strata.names <- as.character(x$stratas[[1]])
      
      for (k in seq_along(params)) {
        
        draw_contour(
          param = params[k],
          tau = taus[k],
          strata.name = strata.names[k],
          derivative = FALSE
        )
        
        draw_contour(
          param = params[k],
          tau = taus[k],
          strata.name = strata.names[k],
          derivative = TRUE
        )
      }
    }
  }
  invisible(NULL)
}
