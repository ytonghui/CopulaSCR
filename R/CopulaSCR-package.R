#' CopulaSCR: Algorithms for fitting copula-based semi-competing risks models
#'
#' The package offers a comprehensive collection of practical and easy-to-use tools for analyzing
#' semi-competing risks with single or multiple intermediate event times and a (possibly) correlated terminal event.
#' The joint distribution of these event times are built upon an Archimedean copula structure.
#' The implemented estimating procedure does not require any parametric assumption on the marginal distributions.
#' The package also includes tools for visualization of semi-competing risks data and
#' simulation from the models.
#' @aliases CopulaSCR-package
#'
#' @references Yu, Tonghui, and Xiang, Liming (2024). Association analysis of multiple intermediate events and dynamic terminal prediction.
#' \emph{Working paper.}
## #' \emph{arXiv}, \bold{1}(0): 0--0.
#' @references Yu, Tonghui, Xiang, Liming, Chen, Chixiang, and Chiou, Sy Han (2024). CopulaSCR:  .
#' \emph{Working paper.}
###' \emph{Journal of XXX}, \bold{105}(5): 1--34.
#'
#' @export
## #' @import copula
## #' @import acopula
## #' @importFrom stats loess model.frame optimize runif
## #' @importFrom prodlim prodlim predictSurvIndividual
#' @importFrom survival survfit Surv
## #' @importFrom quantreg rearrange
#'
## #' @useDynLib CopulaSCR
## #' @useDynLib CopulaSCR, .registration = TRUE
#' @importFrom graphics par text
#' @importFrom stats approx as.formula formula knots mad median update
#' @importFrom stats na.omit predict qnorm quantile rbinom rexp sd stepfun terms
#' @importFrom utils combn
#' @docType package
"_PACKAGE"
NULL


