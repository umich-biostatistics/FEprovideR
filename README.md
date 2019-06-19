
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FEprovideR

## Fixed effects logistic model with high-dimensional parameters

<!-- badges: start -->

<!-- badges: end -->

A stuctured profile likelihood algorithm for the logistic fixed effects
model and an approximate expectation maximization (EM) algorithm for the
logistic mixed effects model.

## Installation

You can install the released version of FEprovideR from
[Github](https://github.com/umich-biostatistics/FEprovideR)
with:

``` r
install.packages("devtools") # you need devtools to install packages from Github
devtools::install_github("umich-biostatistics/FEprovideR")
```

You can install directly from CRAN with:

``` r
install.packages("FEprovideR")
```

## Example

This tutorial simulates a data set to demonstrate the functions provided
by FRprovideR.

``` r
# load the package
library(FEprovideR)

# other imports
library(Matrix)
library(poibin)
library(ggplot2)
```

To simulate a data set, use the following code chunk:

``` r
# Simulate a data set
m <- 500
prov.size <- pmax(round(rnorm(m, 50, 15)),11)
gamma <- rnorm(m, log(3/7), 0.4)
beta <- c(1,0.5,-1)
Y.char <- 'Y'
prov.char <- 'prov.ID'
Z.char <- paste0('z', 1:length(beta))
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
data <- sim.fe.prov(m, prov.size, gamma, beta, Y.char, Z.char, prov.char)
```

This data is also available in the included data sets that come with the
package. To use the included data, run:

``` r
data(hospital)            # raw data
data(hospital_prepared)   # processed data
```

Now, set relevant parameters and fit a model to the prepared
data:

``` r
# a small positive number specifying stopping criterion of Newton-Raphson algorithm
tol <- 1e-5  
# Name input variables and other parameters
Y.char <- 'Y'
prov.char <- 'prov.ID'
Z.char <- paste0('z', 1:3)
data(hospital_prepared) # build in data set
fe.ls <- fe.prov(hospital_prepared, Y.char, Z.char, prov.char, tol) # model fitting
```

Conduct hypothesis tests on the estimated standardized readmission
ratios (SSRs):

``` r
# hypothesis testing
null <- "median"
n <- 10000
alpha <- 0.05
score.fe <- test.fe.prov(hospital_prepared, fe.ls, Y.char, Z.char, prov.char, test="score", null, alpha)
exact.pb <- test.fe.prov(hospital_prepared, fe.ls, Y.char, Z.char, prov.char, test="exact.poisbinom", null, alpha)
exact.bs <- test.fe.prov(hospital_prepared, fe.ls, Y.char, Z.char, prov.char, test="exact.bootstrap", null, alpha, n)
exact.binom <- test.fe.prov(hospital_prepared, fe.ls, Y.char, Z.char, prov.char, test="exact.binom", null="median", alpha)
```

Compute confidence intervals for the estimated SSRs:

``` r
# confidence intervals
confint.df <- confint.fe.prov(fe.ls, parm = "all", level = 0.95, hospital_prepared, Y.char, Z.char, prov.char)
confint.df <- confint.fe.prov(fe.ls, parm = "all", level = 0.90, hospital_prepared, Y.char, Z.char, prov.char)
confint.df <- confint.fe.prov(fe.ls, level = 0.90, data = hospital_prepared, Y.char = Y.char, Z.char = Z.char, prov.char = prov.char)
```

## Funnel plots for SRRs (Standardized readmission ratios)

``` r
# format input data for funnel plot
input.dis <- data.frame(ID=hospital_prepared[hospital_prepared$included==1, prov.char],
                        prob=fe.ls$Exp)
input.prov <- data.frame(SRR=fe.ls$df.prov$SRR, flag=score.fe$flag)
```

Score test based funnel plot:

``` r
target <- c(1)
alphas <- c(0.1, 0.5, 0.01)
input.prov <- data.frame(SRR=fe.ls$df.prov$SRR, flag=score.fe$flag)
funnel.SRR(input.dis, input.prov, target, alphas, type="FE.score")
```
