rm(list = ls())

# test data
n = 500
S = matrix(rnorm(n*n, 0, 1), n, n)


start.time <- Sys.time()
edepths = ED(S)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
