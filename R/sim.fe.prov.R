

sim.fe.prov <- function(m, prov.size, gamma, beta, Y.char, Z.char, prov.char) {
  N <- sum(prov.size) # total number of discharges
  gamma.dis <- rep(gamma, times=prov.size)
  prov <- rep(1:m, times=prov.size) # provider IDs
  Z <- matrix(rnorm(N*length(beta)), ncol=length(beta))
  Y <- rbinom(N, 1, plogis(gamma.dis+Z%*%beta))
  data <- as.data.frame(cbind(Y, prov, Z))
  colnames(data) <- c(Y.char, prov.char, Z.char)
  return(data)
}

