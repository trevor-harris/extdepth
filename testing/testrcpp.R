
rm(list = ls())
gc()

library(microbenchmark)
# Rcpp::sourceCpp("src/test.cpp")

n = 100*100
S = matrix(rnorm(n*100), n, 100)
g = rnorm(n)

microbenchmark(extdepth(S[,1], S), times = 10, unit = "s")

microbenchmark(extdepth_all(S), times = 10, unit = "s")
