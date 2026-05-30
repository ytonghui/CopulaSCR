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
#' @references Yu, Tonghui, and Xiang, Liming. Learning association from multiple intermediate events for dynamic prediction of survival: 
#' an application to cardiovascular disease prognosis,Biometrics, 2026, 82(2): ujag087.
#' 
#' @references Yu, Tonghui, Zhang, Binhui, Xiang, Liming, Ma, Jianwei, and Chen, Chixiang. CopulaSCR: An R Package for the Analysis of Semi-competing Risks Data Using Copula-based Models.
#' \emph{Working paper.}
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
#' @importFrom foreach foreach %dopar%
#' @importFrom stats approxfun optimize
#' @importFrom graphics par text
#' @importFrom stats approx as.formula formula knots mad median update
#' @importFrom stats na.omit predict qnorm quantile rbinom rexp sd stepfun terms
#' @importFrom utils combn
#' @docType package
#' @keywords internal
"_PACKAGE"
NULL
#' @examples
#' \dontrun{
#' if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#' BiocManager::install("graph")
#' BiocManager::install("Rgraphviz")
#' }
#'
#'


