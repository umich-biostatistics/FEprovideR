
#' Simulated readmissions data for 500 hospitals
#'
#' A data set containing simulated readmissions data for 500 hospitals with
#' three continuous covariates. This data needs to be processed with \code{fe.data.prep}.
#'
#' @format A \code{data.frame} with 24438 rows and 5 variables (columns):
#' \describe{
#'   \item{Y}{Indicator for readmission; 1=Yes, 0=No; numeric}
#'   \item{prov.ID}{Provider ID; numeric}
#'   \item{z1}{Simulated covariate 1, numeric}
#'   \item{z2}{Simulated covariate 2, numeric}
#'   \item{z3}{Simulated covariate 3, numeric}
#' }
"hospital"


#' Prepared version of simulated readmissions data for 500 hospitals
#'
#' A data set containing simulated and processed readmissions data for 500 hospitals with
#' three continuous covariates. This is the form of the data needed to use \code{fe.prov}.
#'
#' @format A \code{data.frame} with 24438 rows and 8 variables (columns):
#' \describe{
#'   \item{Y}{Indicator for readmission; 1=Yes, 0=No; numeric}
#'   \item{prov.ID}{Provider ID; numeric}
#'   \item{z1}{Simulated covariate 1, numeric}
#'   \item{z2}{Simulated covariate 2, numeric}
#'   \item{z3}{Simulated covariate 3, numeric}
#'   \item{included}{variable 'included' as an indicator}
#'   \item{no.readm}{providers with no readmission within 30 days}
#'   \item{all.readm}{providers with all readmissions within 30 days}
#' }
"hospital_prepared"
