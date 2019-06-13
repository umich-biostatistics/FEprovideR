#' Prepares data for model fitting (fe.prov)
#'
#' \code{fe.data.prep} prepares the data for model fitting with \code{fe.prov} by
#' taking the data with missing values imputed.
#'
#' @param data a \code{data.frame} including response, provider ID, and covariates, with missing values imputed
#' @param Y.char name of the response variable from \code{data} as a character string
#' @param Z.char names of covariates from \code{data} as vector of character strings
#' @param prov.char name of provider IDs variable as a character string
#' @param cutoff cutoff of provider size as an interger, default value is 10
#'
#' @seealso \code{\link{fe.data.prep}},  \code{\link{fe.prov}},   \code{\link{test.fe.prov}},
#' \code{\link{funnel.SRR}},   \code{\link{confint.fe.prov}}
#'
#' @references He, K., Kalbfleisch, J.D., Li, Y. and Li, Y., 2013. Evaluating hospital
#' readmission rates in dialysis facilities; adjusting for hospital effects. Lifetime data
#' analysis, 19(4), pp.490-512.
#'
#'
#' @examples
#' data(hospital) # build in data set
#' # Name input variables and other parameters
#' cutoff <- 10              # an integer as cutoff of facility (or provider) size with 10 as default
#' alpha <- 0.05             # significance level
#' Y.char <- 'Y'
#' prov.char <- 'prov.ID'
#' Z.char <- paste0('z', 1:3)
#'
#' hospital_prepared <- fe.data.prep(hospital, Y.char, Z.char, prov.char, cutoff) # data preparation
#'
#' @export fe.data.prep
#'
#' @importFrom Matrix rankMatrix

fe.data.prep <- function(data, Y.char, Z.char, prov.char, cutoff=10) {
  #       data: a data frame including response, provider ID, and
  #             covariates, with missing values imputed
  #     Y.char: a character string as name of response variable
  #     Z.char: a vector of character strings as names of covariates
  #  prov.char: a character string as name of variable consisting of provider IDs
  #     cutoff: an integer as cutoff of provider size with 10 as default

  ## check absence of variables
  message("Checking absence of variables ... ")
  Y.ind <- match(Y.char, names(data))
  if (is.na(Y.ind)) {
    stop(paste("Response variable '", Y.char, "' NOT found!", sep=""),call.=F)
  }
  Z.ind <- match(Z.char, names(data))
  if (sum(is.na(Z.ind)) > 0) {
    stop(paste("Covariate(s) '", paste(Z.char[is.na(Z.ind)], collapse="', '"), "' NOT found!", sep=""),call.=F)
  }
  prov.ind <- match(prov.char, names(data))
  if (is.na(prov.ind)) {
    stop(paste("Provider ID '", prov.char, "' NOT found!", sep=""),call.=F)
  }
  message("Checking absence of variables completed!")

  ## check missingness of variables
  message("Checking missingness of variables ... ")
  if (sum(is.na(data[,Y.ind])) > 0) {
    warning(sum(is.na(data[,Y.ind]))," out of ",NROW(data[,Y.ind])," in '",Y.char,"' missing!",immediate.=T,call.=F)
  }
  for (i in 1:length(Z.ind)) {
    if (sum(is.na(data[,Z.ind[i]])) > 0) {
      warning(sum(is.na(data[,Z.ind[i]]))," out of ",NROW(data[,Z.ind[i]])," in '",Z.char[i],"' missing!",immediate.=T,call.=F)
    }
  }
  if (sum(is.na(data[,prov.ind])) > 0) {
    warning(sum(is.na(data[,prov.ind]))," out of ",NROW(data[,prov.ind])," in '",prov.char,"' missing!",immediate.=T,call.=F)
  }
  if (sum(complete.cases(data[,c(Y.ind,Z.ind,prov.ind)]))==NROW(data)) {
    message("Missing values NOT found. Checking missingness of variables completed!")
  } else {
    missingness <- (1 - sum(complete.cases(data[,c(Y.ind,Z.ind,prov.ind)])) / NROW(data)) * 100
    stop(paste(round(missingness,2), "% of all observations are missing!",sep=""),call.=F)
  }
  ## check variation in covariates
  message("Checking variation in covariates ... ")
  ind.novariation <- which(apply(data[, Z.char], 2, var)==0)
  if(length(ind.novariation) > 0){
    stop(paste0("Covariate(s) '", paste(Z.char[ind.novariation], collapse="', '"), "' is/are discharge-invariant!"))
  }
  message("Checking variation in covariates completed!")
  ## check singularity of design matrix
  message("Checking sigularity of design matrix ... ")
  for (i in 2:length(Z.char)) {
    if (rankMatrix(as.matrix(data[,Z.char[1:i]]))[1] < i)
      stop(paste0("Covariate '", Z.char[i], "' is perfectly collinear with the rest!",sep=""),call.=F)
  }
  message("NO multicollinearity in design matrix. Checking sigularity of design matrix completed!")

  data <- data[order(factor(data[,prov.ind])),] # sort data by provider ID
  prov.size <- as.integer(table(data[,prov.ind])) # provider sizes
  prov.size.long <- rep(prov.size,prov.size) # provider sizes assigned to patients
  data$included <- 1 * (prov.size.long > cutoff) # create variable 'included' as an indicator
  warning(sum(prov.size<=cutoff)," out of ",length(prov.size),
          " providers considered small and filtered out!",immediate.=T,call.=F)
  prov.list <- unique(data[data$included==1,prov.ind])   # a reduced list of provider IDs
  prov.no.readm <-      # providers with no readmission within 30 days
    prov.list[sapply(split(data[data$included==1,Y.ind], factor(data[data$included==1,prov.ind])),sum)==0]
  data$no.readm <- 0
  data$no.readm[data[,prov.ind]%in%c(prov.no.readm)] <- 1
  message(paste(length(prov.no.readm),"out of",length(prov.list),
                "remaining providers with no readmission within 30 days."))
  prov.all.readm <-     # providers with all readmissions within 30 days
    prov.list[sapply(split(1-data[data$included==1,Y.ind],factor(data[data$included==1,prov.ind])),sum)==0]
  data$all.readm <- 0
  data$all.readm[data[,prov.ind]%in%c(prov.all.readm)] <- 1
  message(paste(length(prov.all.readm),"out of",length(prov.list),
                "remaining providers with all readmissions within 30 days."))
  message(paste0("After screening, ", round(sum(data[data$included==1,Y.ind])/length(data[data$included==1,Y.ind])*100,2),
                 "% of all discharges were readmitted within 30 days."))
  return(data)
  #       data: a data frame sorted by provider IDs with additional variables 'included', 'no.readm', 'all.readm'
  #             and missing values imputed

}
