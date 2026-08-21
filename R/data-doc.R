#' Example Semi-Competing Risks Dataset with a Single Intermediate Event
#'
#' A simulated dataset for semi-competing risks analysis with a single
#' intermediate event and one terminal event.
#'
#' @format A data frame with 100 rows and 4 variables:
#' \describe{
#'   \item{T1}{Observed time for the intermediate event.}
#'   \item{T2}{Observed time for the terminal event.}
#'   \item{event1}{Event indicator for the intermediate event.}
#'   \item{event2}{Event indicator for the terminal event.}
#' }
#' @usage data(SCRdata)
"SCRdata"

#' Example Semi-Competing Risks Dataset by Treatment Group
#'
#' A simulated dataset for semi-competing risks analysis with a single
#' intermediate event, one terminal event, and a treatment variable.
#'
#' @format A data frame with 500 rows and 5 variables:
#' \describe{
#'   \item{T1}{Observed time for the intermediate event.}
#'   \item{T2}{Observed time for the terminal event.}
#'   \item{event1}{Event indicator for the intermediate event.}
#'   \item{event2}{Event indicator for the terminal event.}
#'   \item{tr}{Treatment group indicator.}
#' }
#' @usage data(SCRdata_by_tr)
"SCRdata_by_tr"

#' Example Semi-Competing Risks Dataset with Multiple Intermediate Events
#'
#' A simulated dataset for semi-competing risks analysis with three
#' intermediate events and one terminal event.
#'
#' @format An object of class \code{mscrData}, implemented as a list with
#' five components:
#' \describe{
#'   \item{T1}{A numeric matrix with 100 rows and 3 columns, containing observed times for three intermediate events.}
#'   \item{T2}{A numeric vector of length 100, containing observed terminal event times.}
#'   \item{event1}{A numeric matrix with 100 rows and 3 columns, containing event indicators for the three intermediate events.}
#'   \item{event2}{A numeric vector of length 100, containing terminal event indicators.}
#'   \item{call}{The function call used to generate the data.}
#' }
#' @usage data(mSCRdata)
"mSCRdata"

#' MIMIC-IV Demo Dataset
#'
#' A demo dataset consisting of 100 patients, derived from a subset of the
#' MIMIC-IV data.
#'
#' @format A data frame with 100 rows and 7 variables:
#' \describe{
#'   \item{person_id}{Patient identifier.}
#'   \item{Digestive}{Recorded time for the digestive event.}
#'   \item{Infection}{Recorded time for the infection event.}
#'   \item{Renal}{Recorded time for the renal event.}
#'   \item{Respiratory}{Recorded time for the respiratory event.}
#'   \item{T_death}{Recorded time of death.}
#'   \item{last_followup}{Recorded time of the last follow-up.}
#' }
#'
#' @source
#' Kallfelz, M. et al. (2021). MIMIC-IV demo data in the OMOP Common Data
#' Model. PhysioNet. RRID: SCR_007345.
#' \doi{10.13026/p1f5-7x35}
#'
#' @usage data(mimiv_demo)
"mimiv_demo"