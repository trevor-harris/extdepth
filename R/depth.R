# Depth Functions

depth <- function(f, fmat) {
  # Computes the depth values of a function with respect to a set of functions (fmat)

  fn = ncol(fmat)
  depth = rep(0, length(f))

  for (row in 1:nrow(fmat)) {
    diff = abs(sum(sign(f[row] - fmat[row,])))
    depth[row] = 1 - (diff / fn)
  }
  return(depth)
}


depth_CDF <- function(f, fmat) {
  # Computes the depth CDF of a function with respect to a set of functions (fmat)

  depths = depth(f, fmat)
  f_obs = nrow(fmat)
  points = seq(1, f_obs) / f_obs

  cdf = sapply(points, function(x) sum(depths <= x) / f_obs)
  return(cdf)
}


point_ED <- function(f1_cdf, f2_cdf) {
  # Calculate if function 1 or function 2 is more extreme
  # returns 1 if f1 is more extreme, -1 if f2 is more extreme, and 0 if equivalent

  for (x in 1:length(f1_cdf)) {
    diff = f1_cdf[x] - f2_cdf[x]
    if (diff != 0.0) {
      return(sign(diff))
    }
  }

  return(0)
}


ED <- function(fmat) {
  # Takes a matrix of functions (each column is a function) and returns there ED ordering
  # from deepest to shallowest

  # find the depth CDF for each function
  dCDF = matrix(0, nrow(fmat), ncol(fmat))
  for (col in 1:ncol(fmat)) {
    dCDF[,col] = depth_CDF(fmat[,col], fmat)
  }

  # for each depth CDF (f1) count the number of depth CDFs (f2) that it's more extreme than
  EDepth = integer(0)
  for (f1 in 1:ncol(dCDF)) {
    gt = 0
    for (f2 in 1:ncol(dCDF)) {
      gt = gt + (point_ED(dCDF[,f1], dCDF[,f2]) > 0)
    }
    EDepth = c(EDepth, gt)
  }

  return(EDepth / nrow(fmat))
}


