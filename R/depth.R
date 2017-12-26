rm(list = ls())

# test data
S = matrix(rnorm(200, 0, 1), 10, 20)
f = rnorm(10, 1, 1)


# pointwise depth (one observation)
obs_depth <- function(obs, set) {
  diff = abs(sum(as.numeric(obs > set) - as.numeric(obs <= set)))
  
  return(1 - (diff/length(set)))
}

# pointwise depth (whole function)
point_depth <- function(func, set) {
  out = rep(0,length(func))
  
  for (row in 1:nrow(set)) {
    out[row] = obs_depth(func, set[row,])
  }
  
  return(out)
}

D = point_depth(f, S)

# depth CDF
depth_CDF <- function(depths) {
  depths.count = length(depths)
  depths.set = seq(1, depths.count) / depths.count
  
  out = sapply(depths.set, function(x) sum(depths <= x) / depths.count)
  return(out)
}

# compute the depth CDF for a set of functions (make into a function)
pdepth = matrix(0, 10, 20)
for (col in 1:ncol(S)) {
  pdepth[,col] = point_depth(S[,col], S)
  pdepth[,col] = depth_CDF(pdepth[,col])
}

Dcdf = depth_CDF(point_depth(f, S))
plot(Dcdf)
