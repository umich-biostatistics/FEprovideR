
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
