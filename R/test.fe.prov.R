#' Hypothesis tests for fe.prov model object
#'
#' \code{test.fe.prov} Conducts hypothesis tests for model parameter estimates.
#' First fit a \code{fe.prov} model object. Go to
#' \href{https://github.com/umich-biostatistics/FEprovideR}{Github} for a tutorial.
#'
#' @param data prepared \code{data.frame}. Use \code{\link{fe.data.prep}}
#' @param fe.ls fitted model object (fit using \code{fe.prov})
#' @param Y.char Y.char name of the response variable from \code{data} as a character string
#' @param Z.char Z.char names of covariates from \code{data} as vector of character strings
#' @param prov.char name of provider IDs variable as a character string
#' @param test string denoting hypothesis test to be conducted. Currently, options
#' include "exact.binom", "exact.poisbinom", "exact.bootstrap", "score". The default
#' is \code{test="score"}
#' @param null use median for null comparison
#' @param alpha alpha level for the CIs
#' @param n number of bootstrap draws
#'
#' @return Returns a \code{data.frame} of the results of the test for each provider
#' with attributes:
#' \itemize{
#'   \item{flag:} Either "1" for p<alpha/2, "0" p<=1-alpha/2 and p<alpha/2, or "-1" for neither
#'   \item{p:} p-value for the hypothesis test of the model parameter
#' }
#'
#' @seealso \code{\link{fe.data.prep}},  \code{\link{fe.prov}},
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
#' # Hypothesis tests
#' null = "median"
#' alpha = 0.05
#' score.fe <- test.fe.prov(hospital_prepared, fe.ls, Y.char, Z.char,
#'                          prov.char, test="score", null, alpha)
#'
#' @export test.fe.prov
#'
#' @importFrom poibin ppoibin
#' @import stats


test.fe.prov <- function(data, fe.ls, Y.char, Z.char, prov.char, test="score", null="median", alpha=0.05, n=10000) {

  if (!(test %in% c("exact.binom", "exact.poisbinom", "exact.bootstrap", "score")))
    stop("Argument 'test' NOT as required!",call.=F)
  data <- data[data$included==1, c(Y.char, Z.char, prov.char)]
  gamma <- fe.ls$df.prov$gamma; beta <- fe.ls$beta
  gamma.null <- ifelse(null=="median", median(gamma),
                       ifelse(class(null)=="numeric", null[1],
                              stop("Argument 'null' NOT as required!",call.=F)))
  if (test=="exact.bootstrap") {
    exact.bootstrap <- function(df, n) {
      probs <- plogis(gamma.null + unname(as.matrix(df[, Z.char])) %*% beta)
      obs <- sum(df[,Y.char])
      sums <- colSums(matrix(rbinom(n=length(probs)*n, size=1, prob=rep(probs,times=n)), ncol=n))
      p <- (sum(sums>obs)+ 0.5*sum(sums==obs))/n
      flag <- ifelse(p<alpha/2, 1, ifelse(p<=1-alpha/2, 0, -1))
      p.val <- 2 * min(p, 1-p)
      return(c(flag, p.val))
    }
    results <- sapply(by(data, data[,prov.char],identity),
                      FUN=function(x) exact.bootstrap(x, n))
    return(data.frame(flag=factor(results[1,]), p=results[2,]))
  }  else if (test=="score") {
    data$probs <- plogis(gamma.null + unname(as.matrix(data[, Z.char])) %*% beta)
    z.score <- sapply(split(data[,Y.char]-data$probs,data[,prov.char]),sum) /
      sqrt(sapply(split(data$probs*(1-data$probs),data[,prov.char]),sum))
    p <- pnorm(z.score, lower.tail=F)
    flag <- ifelse(p<alpha/2, 1, ifelse(p<=1-alpha/2, 0, -1))
    p.val <- 2 * pmin(p, 1-p)
    return(data.frame(flag=factor(flag), p=p.val, row.names=unique(data[, prov.char])))
  } else if (test=="exact.poisbinom") {
    exact.poisbinom <- function(df) {
      probs <- plogis(gamma.null + unname(as.matrix(df[, Z.char])) %*% beta)
      obs <- sum(df[,Y.char])
      p <- 1 - ppoibin(obs, probs) + 0.5*dpoibin(obs, probs)
      flag <- ifelse(p<alpha/2, 1, ifelse(p<=1-alpha/2, 0, -1))
      p.val <- 2 * min(p, 1-p)
      return(c(flag, p.val))
    }
    results <- sapply(by(data, data[,prov.char],identity),
                      FUN=function(x) exact.poisbinom(x))
    return(data.frame(flag=factor(results[1,]), p=results[2,]))
  } else if (test=="exact.binom") {
    exact.binom <- function(df) {
      probs <- plogis(gamma.null + unname(as.matrix(df[, Z.char])) %*% beta)
      obs <- sum(df[,Y.char])
      p <- 1 - pbinom(obs, size=length(probs), prob=mean(probs)) +
        0.5*dbinom(obs, size=length(probs), prob=mean(probs))
      flag <- ifelse(p<alpha/2, 1, ifelse(p<=1-alpha/2, 0, -1))
      p.val <- 2 * min(p, 1-p)
      return(c(flag, p.val))
    }
    results <- sapply(by(data, data[,prov.char],identity),
                      FUN=function(x) exact.binom(x))
    return(data.frame(flag=factor(results[1,]), p=results[2,]))
  }
}
