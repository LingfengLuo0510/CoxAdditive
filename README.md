# CoxAdditive
# The example code to install this package and run the code:
library(Rcpp)

library(RcppArmadillo)

library(mvtnorm)

library(Matrix)

library(MASS)

library(splines)

library(devtools)

install_github("https://github.com/LingfengLuo0510/CoxAdditive")

s       <- simul(p=100,N=500,knot = 6)  

x       <- CoxMain(s$delta,s$z,s$time,s$p,s$N, s$knot,maxit = 300)

y       <- CoorDesc(s$delta,s$z,s$time,s$p,s$N, s$knot, maxit = 300)

x$time

y$time

FPFNSeSpLik(s$TrueBeta,x$beta)

FPFNSeSpLik(s$TrueBeta,y$beta)

