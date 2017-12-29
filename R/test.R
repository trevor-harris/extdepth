rm(list = ls())

source('R/depth.R')
source('R/central_regions.R')

# timing
timeit <- function(func) {
  start.time <- Sys.time()
  func
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time.taken
}

# test data
n = 100
S = matrix(rnorm(n*n), n, n)

edepths = ED(S)
cr = central_region(ED, 0.05)

plot(S[cr$median,], type = 'l', col="blue")
lines(S[cr$lower,], col="red")
lines(S[cr$upper,], col="red")

# timeit(ED(S))

