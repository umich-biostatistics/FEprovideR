#' Compute confidence intervals for fitted model
#'
#' \code{confint.fe.prov} computes the (1-alpha)\% confidence intervals for the fixed
#' effect parameter estimates. Go to
#' \href{https://github.com/umich-biostatistics/FEprovideR}{Github} for a tutorial.
#'
#' @param object fitted model object (fit using \code{fe.prov})
#' @param parm parameter names. Since their are so many parameters, the default is
#'  \code{"all"}, which is currently the only option available
#' @param level confidence level (default is \code{0.95})
#' @param data prepared \code{data.frame}. Use \code{\link{fe.data.prep}} to prepare the raw data
#' @param Y.char Y.char name of the response variable from \code{data} as a character string
#' @param Z.char Z.char names of covariates from \code{data} as vector of character strings
#' @param prov.char name of provider IDs variable as a character string
#' @param ... extra arguments to be passed to confint
#'
#' @return Returns a \code{data.frame} of gamma and SRR lower and upper CI bounds. Each row is a
#' parameter, each column gives a different bound.
#'
#' @seealso \code{\link{fe.data.prep}},  \code{\link{fe.prov}},   \code{\link{test.fe.prov}},
#' \code{\link{funnel.SRR}},   \code{\link{confint.fe.prov}}
#'
#' @references He, K., Kalbfleisch, J.D., Li, Y. and Li, Y., 2013. Evaluating hospital
#' readmission rates in dialysis facilities; adjusting for hospital effects. Lifetime data
#' analysis, 19(4), pp.490-512.
#'
#' @examples
#' # Name input variables and other parameters
#' # a small positive number specifying stopping
#' # criterion of Newton-Raphson algorithm
#' tol <- 1e-5
#' Y.char <- 'Y'
#' prov.char <- 'prov.ID'
#' Z.char <- paste0('z', 1:3)
#' data(hospital_prepared) # build in data set
#' fe.ls <- fe.prov(hospital_prepared, Y.char, Z.char, prov.char, tol) # model fitting
#'
#' # confidence intervals
#' confint.fe.prov(fe.ls, parm = "all", level = 0.95, hospital_prepared, Y.char, Z.char, prov.char)
#'
#' @export confint.fe.prov
#' @importFrom poibin ppoibin
#' @importFrom poibin dpoibin
#' @export

confint.fe.prov <- function(object, parm = "all", level = 0.95, data, Y.char, Z.char, prov.char,...) {
  if(!(parm == "all")) stop("This method is only implemented for all parameters.
                            You cannot get CI's for a subset.")
  fe.ls <- object
  alpha <- 1 - level
  data <- data[data$included==1, ]
  df.prov <- fe.ls$df.prov
  gamma <- df.prov$gamma; names(gamma) <- rownames(df.prov)
  beta <- fe.ls$beta
  max.gamma <- norm(as.matrix(gamma[is.finite(gamma)]),"I")

  CL.finite <- function(df) {
    UL.gamma <- function(Gamma)
      ppoibin(Obs-1,plogis(Gamma+Z.beta))+0.5*dpoibin(Obs-1,plogis(Gamma+Z.beta))-alpha/2
    LL.gamma <- function(Gamma)
      1-ppoibin(Obs,plogis(Gamma+Z.beta))+0.5*dpoibin(Obs-1,plogis(Gamma+Z.beta))-alpha/2
    prov <- ifelse(length(unique(df[,prov.char]))==1, unique(df[,prov.char]),
                   stop("Number of providers involved NOT equal to one!"))
    Z.beta <- as.matrix(df[,Z.char])%*%beta
    Obs <- df.prov[prov, "Obs"]; Exp <- df.prov[prov, "Exp"]
    gamma.lower <- uniroot(LL.gamma, gamma[prov]+c(-5,0))$root
    gamma.upper <- uniroot(UL.gamma, gamma[prov]+c(0,5))$root
    SRR.lower <- sum(plogis(gamma.lower+Z.beta)) / Exp
    SRR.upper <- sum(plogis(gamma.upper+Z.beta)) / Exp
    return(c(gamma.lower, gamma.upper, SRR.lower, SRR.upper))
  }
  CL.no.readm <- function(df) {
    prov <- ifelse(length(unique(df[,prov.char]))==1, unique(df[,prov.char]),
                   stop("Number of providers involved NOT equal to one!"))
    Z.beta <- as.matrix(df[,Z.char])%*%beta
    max.Z.beta <- norm(Z.beta, "I")
    gamma.upper <- uniroot(function(x) prod(plogis(-x-Z.beta))/2-alpha,
                           (10+max.Z.beta)*c(-1,1))$root
    SRR.upper <- sum(plogis(gamma.upper+Z.beta)) / df.prov[prov, "Exp"]
    return(c(-Inf, gamma.upper, 0, SRR.upper))
  }
  CL.all.readm <- function(df) {
    prov <- ifelse(length(unique(df[,prov.char]))==1, unique(df[,prov.char]),
                   stop("Number of providers involved NOT equal to one!"))
    Z.beta <- as.matrix(df[,Z.char])%*%beta
    max.Z.beta <- norm(Z.beta, "I")
    Exp <- df.prov[prov, "Exp"]; SRR <- df.prov[prov, "SRR"]
    gamma.lower <- uniroot(function(x) prod(plogis(x+Z.beta))/2-alpha,
                           (10+max.Z.beta)*c(-1,1))$root
    SRR.lower <- sum(plogis(gamma.lower+Z.beta)) / Exp
    return(c(gamma.lower, Inf, SRR.lower, SRR))
  }
  confint.finite <- sapply(by(data[(data$no.readm==0) & (data$all.readm==0),], data[(data$no.readm==0) & (data$all.readm==0),prov.char],identity),
                           FUN=function(df) CL.finite(df))
  confint.no.readm <- sapply(by(data[data$no.readm==1,], data[data$no.readm==1,prov.char],identity),
                             FUN=function(df) CL.no.readm(df))
  confint.all.readm <- sapply(by(data[data$all.readm==1,], data[data$all.readm==1,prov.char],identity),
                              FUN=function(df) CL.all.readm(df))
  confint.df <- as.data.frame(t(cbind(confint.finite, confint.no.readm, confint.all.readm)))
  names(confint.df) <- c("gamma.lower", "gamma.upper", "SRR.lower", "SRR.upper")
  return(confint.df[order(rownames(confint.df)),])
}
