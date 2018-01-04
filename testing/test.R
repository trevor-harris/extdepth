
rm(list = ls())
gc()

library(microbenchmark)
library(extdepth)

# test data
n = 100
S = matrix(rnorm(2*n*n), 2*n, n)

microbenchmark(dCDF(S[,1], S), times = 100)

microbenchmark(depth(S[,1], S), times = 100)


# without compiler
woc = microbenchmark(ED(S), times = 100)

# with compiler
ED = compiler::cmpfun(ED)
wc = microbenchmark(ED(S), times = 100)

woc
wc

edepths = ED(S)
cr = central_region(edepths, 0.05)

plot(S[cr$median,], type = 'l', col="blue")
lines(S[cr$lower,], col="red")
lines(S[cr$upper,], col="red")





