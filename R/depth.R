rm(list = ls())


# pointwise depth (one observation)
obs_depth <- function(obs, vec) {
  diff = abs(sum(as.numeric(obs > vec) - as.numeric(obs <= vec)))
  depth = 1 - (diff / length(vec))
  return(depth)
}

# pointwise depth (whole function)
point_depth <- function(f, fmat) {
  depth = rep(0, length(f))

  for (row in 1:nrow(fmat)) {
    depth[row] = obs_depth(f[row], fmat[row,])
  }
  return(depth)
}

# depth CDF
depth_CDF <- function(depths) {
  len = length(depths)
  points = seq(1, len) / len

  cdf = sapply(points, function(x) sum(depths <= x) / len)
  return(cdf)
}

point_ED <- function(f_cdf, g_cdf) {
  # returns 1 if f is more extreme and -1 if g is more extreme
  for (x in 1:length(f_cdf)) {
    diff = f_cdf[x] - g_cdf[x]
    if (diff != 0.0) {
      return(sign(diff))
    }
  }
  # incase neither is more extreme
  return(0)
}

ED <- function(fmat) {
  # find the depth CDF for each function
  dCDF = matrix(0, nrow(fmat), ncol(fmat))
  for (col in 1:ncol(fmat)) {
    dCDF[,col] = point_depth(fmat[,col], fmat)
    dCDF[,col] = depth_CDF(dCDF[,col])
  }

  # for each depth CDF (f1) count the number of depth CDFs (f2) that it's more extreme than
  gt.vec = integer(0)
  for (f1 in 1:ncol(dCDF)) {
    gt = 0
    for (f2 in 1:ncol(dCDF)) {
      gt = gt + (point_ED(dCDF[,f1], dCDF[,f2]) > 0)
    }
    gt.vec = c(gt.vec, gt)
  }

  return(rev(order(gt.vec)))
}

# test data
S = matrix(rnorm(10000, 0, 1), 1000, 10)
S =
ED(S)



