
rm(list = ls())
gc()

library(microbenchmark)
library(extdepth)

n = 100*100
S = matrix(rnorm(n*100), n, 100)

f1 = S[,1]
f2 = S[,2]

cdf1 = dCDF(depth(f1, S))
cdf2 = dCDF(depth(f2, S))
# 
# dCDF(f1, S)
# depth_CDF(f1, S)

# cdf calc check
f1depth = depth(f1, S)
microbenchmark(dCDF(f1depth), times = 20, unit = "ms")
microbenchmark(depth_CDFr(f1, S), times = 20, unit = "ms")

# pointwise ed check
t = 10000
microbenchmark(ed_compare(cdf1, cdf2), times = t, unit = "us")
microbenchmark(point_EDr(cdf1, cdf2), times = t, unit = "us")

# full check
microbenchmark(ED(S), times = 10, unit = "s")
microbenchmark(EDr(S), times = 10, unit = "s")

