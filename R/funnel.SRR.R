
funnel.SRR <- function(input.dis, input.prov, target=1, alphas=c(0.1, 0.05, 0.01), type="FE.score", file, sigma.b=NULL){
  # input.dis: a data frame consisting of discharge-specific inputs and provider ID
  # input.prov: a data frame consisting of provider-specific inputs ordered by provider ID
  # target
  # alphas
  # type
  # file
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
  ggsave(file, width=8, height=8)
}
