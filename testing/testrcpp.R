
rm(list = ls())
gc()

library(microbenchmark)
Rcpp::sourceCpp("src/test.cpp")

n = 100*1000
S = matrix(rnorm(n*100), n, 100)

microbenchmark(fast_ED(S), times = 1, unit = "s")
microbenchmark(EDr(S), times = 1, unit = "s")

1/n

####
f1 = S[,1]
f2 = S[,2]

cdf1 = dCDF(depth(f1, S))
cdf2 = dCDF(depth(f2, S))

# cdf calc check
f1depth = depth(f1, S)
f1depth[which.min(f1depth)]
microbenchmark(dCDF(f1depth), times = 20, unit = "ms")
microbenchmark(depth_CDFr(f1, S), times = 20, unit = "ms")

# pointwise ed check
microbenchmark(ed_compare(cdf1, cdf2), times = 10000, unit = "us")
microbenchmark(point_EDr(cdf1, cdf2), times = 10000, unit = "us")

# full check
# microbenchmark(ED(S), times = 10, unit = "s")
microbenchmark(fast_ED(S), times = 1, unit = "s")
microbenchmark(EDr(S), times = 1, unit = "s")

f1 = S[,1]
fmat = S[,-1]

ed1 = ED_f(f1,S)
microbenchmark(ED_f(f1,S), times = 10, unit = "s")

