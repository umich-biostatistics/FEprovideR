
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FEprovideR

## Fixed effects logistic model with high-dimensional parameters

<!-- badges: start -->

<!-- badges: end -->

A stuctured profile likelihood algorithm for the logistic fixed effects
model and an aproximate EM algorithm for the logistic mixed effects
model.

## Installation

You can install the released version of FEprovideR from
[Github](https://github.com/umich-biostatistics/FEprovideR)
with:

``` r
install.packages("devtools") # you need devtools to install packages from Github
devtools::install_github("umich-biostatistics/FEprovideR")
```

## Example

This tutorial simulates a data set to demonstrate the functions provided
by FRprovideR.

``` r
# load the package
#library(FEprovideR)
## basic example code

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
data <- sim.fe.prov(m, prov.size, gamma, beta, Y.char, Z.char, prov.char)
```

Now, set relevant parameters and fit a model to the prepared
data:

``` r
cutoff <- 10              # an integer as cutoff of facility (or provider) size with 10 as default
tol <- 1e-5               # a small positive number specifying stopping criterion of Newton-Raphson algorithm
n <- 10000                # resample size, 10000 as default
alpha <- 0.05             # significance level
data.prepared <- fe.data.prep(data, Y.char, Z.char, prov.char, cutoff) # data preparation
fe.ls <- fe.prov(data.prepared, Y.char, Z.char, prov.char, tol) # model fitting
```

Conduct hypothesis tests on the estimated standardized readmission
ratios (SSR):

``` r
# hypothesis testing
null <- "median"
score.fe <- test.fe.prov(data.prepared, fe.ls, Y.char, Z.char, prov.char, test="score", null, alpha)
exact.pb <- test.fe.prov(data.prepared, fe.ls, Y.char, Z.char, prov.char, test="exact.poisbinom", null, alpha)
exact.bs <- test.fe.prov(data.prepared, fe.ls, Y.char, Z.char, prov.char, test="exact.bootstrap", null, alpha, n)
exact.binom <- test.fe.prov(data.prepared, fe.ls, Y.char, Z.char, prov.char, test="exact.binom", null="median", alpha)
```

Compute confidence intervals for the estimated SSR’s:

``` r
# confidence intervals
confint.df <- confint.fe.prov(data.prepared, fe.ls, Y.char, Z.char, prov.char, alpha)
confint.df <- confint.fe.prov(data.prepared, fe.ls, Y.char, Z.char, prov.char, alpha = 0.1)
```

## Funnel plots for SSR

``` r
# Prepare funnel plot inputs
target <- 1
alphas <- c(0.01, 0.05, 0.1)
input.dis <- data.frame(ID=data.prepared[data.prepared$included==1, prov.char], prob=fe.ls$Exp)
```

Score test based funnel plot:

``` r
# file is saved to working directory
input.prov <- data.frame(SRR=fe.ls$df.prov$SRR, flag=score.fe$flag)
funnel.SRR(input.dis, input.prov, target, alphas, type="FE.score", file="./SRR_funnel_fe_score.pdf")
```

Exact test based funnel plot:

``` r
input.prov <- data.frame(SRR=fe.ls$df.prov$SRR, flag=exact.pb$flag)
funnel.SRR(input.dis, input.prov, target, alphas, type="FE.exact", file="./SRR_funnel_fe_exact.pdf")
```
