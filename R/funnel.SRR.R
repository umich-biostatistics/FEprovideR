#'
#' Funnel plot for SRR (standardized readmission ratios)
#'
#' \code{funnel.SRR} produces and returns funnel plots for the analysis using discharge-specific
#' and patient-specific inputs with provider ID.
#'
#' @param input.dis a \code{data.frame} consisting of discharge-specific inputs and provider ID
#' @param input.prov a \code{data.frame} consisting of provider-specific inputs and provider ID
#' @param target target standardized readmission ratio (SRR)
#' @param alphas numeric vector of alpha levels of interest
#' @param type string of length one containing the type of test performed. Currently options
#' include "score", "exact", "FE.score", "FE.exact", "FERE.score", "FERE.exact"
#' @param sigma.b sigma for random effects. Should only have value other than null
#' if prefix "FERE." specified
#' in \code{type=} argument
#'
#' @return Returns a \code{ggplot} object.
#'
#'
#' @seealso \code{\link{fe.data.prep}},  \code{\link{fe.prov}},   \code{\link{test.fe.prov}},
#' \code{\link{confint.fe.prov}}, \code{ggplot2}
#'
#' @references He, K., Kalbfleisch, J.D., Li, Y. and Li, Y., 2013. Evaluating hospital
#' readmission rates in dialysis facilities; adjusting for hospital effects. Lifetime data
#' analysis, 19(4), pp.490-512.
#'
#' @examples
#' # Name input variables and other parameters
#' tol <- 1e-5               # a small positive number specifying stopping criterion of Newton-Raphson algorithm
#' Y.char <- 'Y'
#' prov.char <- 'prov.ID'
#' Z.char <- paste0('z', 1:3)
#' data(hospital_prepared) # build in data set
#' fe.ls <- fe.prov(hospital_prepared, Y.char, Z.char, prov.char, tol) # model fitting
#'
#'
#' # Hypothesis tests
#' null = "median"
#' alpha <- 0.05             # significance level
#' score.fe <- test.fe.prov(hospital_prepared, fe.ls, Y.char, Z.char, prov.char, test="score", null, alpha)
#'
#' # format input data for funnel plot
#' input.dis <- data.frame(ID=hospital_prepared[hospital_prepared$included==1, prov.char], prob=fe.ls$Exp)
#' input.prov <- data.frame(SRR=fe.ls$df.prov$SRR, flag=score.fe$flag)
#'
#' # render funnel plot
#' target <- c(1)
#' alphas = c(0.1, 0.05, 0.01)
#' funnel.SRR(input.dis, input.prov, target, alphas, type="FE.score")
#'

funnel.SRR <- function(input.dis, input.prov, target=1, alphas=c(0.1, 0.05, 0.01), type="FE.score", sigma.b=NULL){
  # input.dis: a data frame consisting of discharge-specific inputs and provider ID
  # input.prov: a data frame consisting of provider-specific inputs ordered by provider ID
  # target
  # alphas
  # type
  # sigma.b
  if (length(unique(input.dis$ID))!=NROW(input.prov))
    stop("Number of unique provider IDs NOT equal to length of indicator vector!",.call=F)
  n <- length(alphas); m <- NROW(input.prov); alphas <- alphas[order(alphas)]
  mean.obs.prov <- as.numeric(sapply(split(input.dis$prob,input.dis$ID),sum))
  var.obs.prov <- as.numeric(sapply(split(input.dis$prob*(1-input.dis$prob),input.dis$ID),sum))
  if (grepl("FE.",type)) {
    se.SRR <- sqrt(var.obs.prov) / mean.obs.prov
  }  else if (grepl("FERE.",type)) {
    if (is.null(sigma.b)) stop("Arugument 'sigma.b' NOT as required!")
    se.SRR <- sqrt(var.obs.prov * (1+sigma.b*var.obs.prov)) / mean.obs.prov
  }
  order.prec <- order(se.SRR, decreasing=T)
  input.prov <- input.prov[order.prec,]; se.SRR <- se.SRR[order.prec]
  data <- data.frame(precision=se.SRR^{-2}, indicator=input.prov$SRR, flag=input.prov$flag)
  ctrl.limits <- data.frame(precision=rep(se.SRR^{-2}, each=n), alpha=rep(alphas, times=m),
                            exp=rep(mean.obs.prov[order.prec], each=n))
  if (grepl("score",type)) {
    ctrl.limits$upper <- target+qnorm(1-ctrl.limits$alpha/2)*sqrt(1/ctrl.limits$precision)
    ctrl.limits$lower <- target-qnorm(1-ctrl.limits$alpha/2)*sqrt(1/ctrl.limits$precision)
  } else if (grepl("exact",type)) {
    CL.obs <- function(df) {
      aux <- function(alpha) {
        # lower CL for obs
        o <- qpoibin(alpha/2, df$prob)
        o <- ifelse(ppoibin(o-1,df$prob)+0.5*dpoibin(o,df$prob)>=alpha/2, o, o+1)
        lambda <- (dpoibin(o,df$prob)+2*ppoibin(o-1,df$prob)-alpha) /
          (dpoibin(o,df$prob)+dpoibin(o-1,df$prob))
        lower <- pmax(o-lambda,0)
        # upper CL for obs
        o <- qpoibin(1-alpha/2, df$prob)
        o <- ifelse(ppoibin(o-1,df$prob)+0.5*dpoibin(o,df$prob)>=1-alpha/2, o, o+1)
        lambda <- (dpoibin(o,df$prob)+2*ppoibin(o-1,df$prob)-2+alpha) /
          (dpoibin(o-1,df$prob)+dpoibin(o,df$prob))
        upper <- pmin(o-lambda, length(df$prob))
        return(c(lower, upper))
      }
      return(as.vector(t(sapply(alphas, FUN=aux))))
    }
    mat <- unname(sapply(by(input.dis,input.dis[,"ID"],identity),FUN=function(df) CL.obs(df)))[,order.prec]
    ctrl.limits$lower <- as.vector(mat[1:n,]) / ctrl.limits$exp
    ctrl.limits$upper <- as.vector(mat[(n+1):(2*n),]) / ctrl.limits$exp

  }
  ctrl.limits$alpha <- factor(ctrl.limits$alpha)
  xmax <- max(se.SRR^{-2}); ymax <- max(input.prov$SRR)
  labs.color <- paste0(levels(input.prov$flag)," (",round((summary(input.prov$flag)/m*100),digits=2),")")
  labs.linetype <- paste0((1-alphas)*100,"%")
  values.linetype <- c('dashed','dotted','dotdash','longdash','twodash')[1:length(alphas)]
  values.linetype[alphas==0.05] <- 'solid'

  ggplot() + theme_classic() +
    theme(legend.justification=c(1,1), legend.position=c(1,1), legend.title=element_text(face="bold")) +
    scale_x_continuous(name="Effective Provider Size", limits=c(0, xmax)) +
    scale_y_continuous(name="SRR", breaks=round(seq(0,ymax,by=1),1), limits=c(0, ymax)) +
    geom_point(data=data, aes(x=precision, y=indicator, color=flag), shape=20, size=0.9, alpha=1) + # "#619CFF"
    scale_color_manual(values=c('darkblue','#E69F00', 'darkred'), labels=labs.color) +
    geom_line(data=ctrl.limits, aes(x=precision, y=lower, group=alpha, linetype=alpha), size=.6) +
    geom_line(data=ctrl.limits, aes(x=precision, y=upper, group=alpha, linetype=alpha), size=.6) +
    scale_linetype_manual(values=values.linetype, labels=labs.linetype) +
    labs(linetype="Ctrl Limits", color="Flagging (%)") +
    guides(color=guide_legend(order=1), linetype=guide_legend(reverse=TRUE, order=2)) +
    geom_hline(yintercept=target, size=.6, linetype="dashed") # color="#F8766D"
}
