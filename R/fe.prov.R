#'
#' Fit logistic fixed-effect model with high-dimensional predictors
#'
#' \code{fe.prov} fits a fixed-effect logistic model using structured profile
#' likelihood algorithm. Standardized readmission ratios (SRRs) are also computed.
#' Go to \href{https://github.com/umich-biostatistics/FEprovideR}{Github} for
#' a tutorial.
#'
#'
#' @param data prepared \code{data.frame}. Use \code{\link{fe.data.prep}} to prepare the raw data
#' @param Y.char name of the response variable from \code{data} as a character string
#' @param Z.char names of covariates from \code{data} as vector of character strings
#' @param prov.char name of provider IDs variable as a character string
#' @param tol tolerance level for convergence. Default is \code{1e-5}
#' @param null use median for null comparison
#'
#' @return An object of class \code{fe.prov}, which is just a \code{List} object with the following named elements:
#' \itemize{
#'   \item \code{beta:} a vector of fixed effect estimates
#'   \item \code{Obs:} a vector of responses for included providers
#'   \item \code{Exp:} a vector of expected probabilities of readmission within 30 days of discharge
#'   \item \code{iter:} number of iterations needed for convergence
#'   \item \code{beta.max.diff:} value of the stopping criterion
#'   \item \code{df.prov:}
#' }
#' \code{df.prov} is a \code{data.frame} of provider-level information with the following items:
#' \itemize{
#'   \item \code{Obs:} provider-level observed number of readmissions within 30 days
#'   \item \code{Exp:} expected number of readmissions within 30 days
#'   \item \code{SRR:} standardized readmission ratios for each hospital
#'   \item \code{gamma:} a vector of provider effect estimates for included hospitals
#' }
#'
#'
#' @seealso \code{\link{fe.data.prep}},   \code{\link{test.fe.prov}},
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
#' @export fe.prov
#'

fe.prov <- function(data, Y.char, Z.char, prov.char, tol=1e-5, null="median"){
  #       data: a data frame sorted by providers with additional variable 'included',
  #             with missing values imputed
  #     Y.char: a character string as name of response variable
  #     Z.char: a vector of character strings as names of covariates
  #  prov.char: a character string as name of variable consisting of provider IDs
  #        tol: a small positive number specifying stopping criterion of Newton-Raphson algorithm

  data <- data[data$included==1,]
  n.prov <- sapply(split(data[, Y.char], data[, prov.char]), length)
  gamma.prov <- rep(0, length(n.prov))
  Z <- as.matrix(data[,Z.char])
  beta <- rep(0, NCOL(Z))

  max.iter <- 10000
  iter <- 0
  bound <- 20
  beta.max.diff <- 100 # initialize stop criterion
  # commented cat for CRAN standards
  # cat("Implementing Newton-Raphson algorithm for fixed provider effects model ...")

  while (iter<=max.iter & beta.max.diff>=tol) {
    iter <- iter + 1
    # cat(paste0("\n Iter ",iter,":")) # remove pringing to console
    # provider effect
    gamma.obs <- rep(gamma.prov, n.prov)
    Z.beta <- Z%*%beta
    p <- c(plogis(gamma.obs+Z.beta))
    q <- 1-p; pq <- p*q
    gamma.prov <- sapply(split(data[,Y.char]-p, data[,prov.char]), sum) /
      sapply(split(pq, data[,prov.char]), sum) + gamma.prov
    gamma.prov <- pmin(pmax(gamma.prov, -bound), bound)
    gamma.obs <- rep(gamma.prov, n.prov)

    # covar parameters
    p <- plogis(gamma.obs+Z.beta)
    q <- 1-p; pq <- p*q
    beta.score <- t(Z)%*%(data[,Y.char]-p)
    beta.info <- t(Z)%*%(c(pq)*Z)
    beta.new <- beta + as.numeric(solve(beta.info)%*%beta.score) # pmin(pmax(, -bound), bound)
    beta.max.diff <- max(abs(beta-beta.new)) # stopping criterion
    beta <- beta.new
    # removed printing to console for CRAN standards
    # cat(paste0(" Inf norm of running diff in est covar coef is ",round(beta.max.diff,digits=8),";"))
  }
  # removed printing to console
  # cat(paste("\n Newton-Raphson algorithm converged after",iter,"iterations! \n"))
  gamma.prov[gamma.prov==bound] <- Inf; gamma.prov[gamma.prov==-bound] <- -Inf
  Z.beta <- as.matrix(data[,Z.char]) %*% beta
  gamma.null <- ifelse(null=="median", median(gamma.prov),
                       ifelse(class(null)=="numeric", null[1],
                              stop("Argument 'null' NOT as required!",call.=F)))
  Exp <- as.numeric(plogis(gamma.null+Z.beta)) # expected prob of readm within 30 days of discharge under null

  df.prov <- data.frame(Obs=sapply(split(data[,Y.char],data[,prov.char]),sum),
                        Exp=sapply(split(Exp,data[,prov.char]),sum))
  df.prov$SRR <- df.prov$Obs / df.prov$Exp
  df.prov$gamma <- gamma.prov
  return(structure(list(beta=beta, Obs=data[, Y.char], Exp=Exp, df.prov=df.prov,
                        iter = iter, beta.max.diff = beta.max.diff), class = "fe.prov"))
  #    beta: a vector of fixed effect estimates
  #     Obs: a vector of responses for included providers
  #     Exp: a vector of expected probs of readmission within 30 days of discharge
  # df.prov: a data frame of provider-level number of observed number of readmissions within 30 days
  #          expected number of readmissions within 30 days, a vector of SRRs, and a vector of
  #          provider effect estimates for included providers (considered as a fixed effect)
}


