#' Compute confidence intervals for fitted model
#'
#' \code{confint.fe.prov} computes the (1-alpha)% confidence interval for the fixed
#' effect parameter estimates
#'
#' @param data prepared \code{data.frame}. Use \code{\link{fe.data.prep}}
#' @param fe.ls fitted model object (fit using \code{fe.prov})
#' @param Y.char Y.char name of the response variable from \code{data} as a character string
#' @param Z.char Z.char names of covariates from \code{data} as vector of character strings
#' @param prov.char name of provider IDs variable as a character string
#' @param alpha alpha level for the CIs
#'
#' @return Returns a \code{data.frame} of gamma and SRR lower and upper CI bounds. Each row is a
#' parameter, each column gives a different bound.
#'
#' @seealso \code{\link{fe.data.prep}},  \code{\link{fe.prov}},   \code{\link{test.fe.prov}},
#' \code{\link{funnel.SRR}},   \code{\link{confint.fe.prov}}
#'
#' @examples
#' ## data simulation
#' m <- 500
#' prov.size <- pmax(round(rnorm(m, 50, 15)),11)
#' gamma <- rnorm(m, log(3/7), 0.4)
#' beta <- c(1,0.5,-1)
#' Y.char <- 'Y'
#' prov.char <- 'prov.ID'
#' Z.char <- paste0('z', 1:length(beta))
#' data <- sim.fe.prov(m, prov.size, gamma, beta, Y.char, Z.char, prov.char)
#'
#' library(Matrix)
#' library(poibin)
#' library(ggplot2)
#' cutoff <- 10              # an integer as cutoff of facility (or provider) size with 10 as default
#' tol <- 1e-5               # a small positive number specifying stopping criterion of Newton-Raphson algorithm
#' n <- 10000                # resample size, 10000 as default
#' alpha <- 0.05             # significance level
#' data.prepared <- fe.data.prep(data, Y.char, Z.char, prov.char, cutoff) # data preparation
#' fe.ls <- fe.prov(data.prepared, Y.char, Z.char, prov.char, tol) # model fitting
#'
#' # confidence intervals
#' confint.fe.prov(data.prepared, fe.ls, Y.char, Z.char, prov.char, alpha)
#'
confint.fe.prov <- function(data, fe.ls, Y.char, Z.char, prov.char, alpha) {

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
