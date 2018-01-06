# Depth Functions

#' @export
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

#' @export
depth_CDF <- function(f, fmat) {
  # Computes the depth CDF of a function with respect to a set of functions (fmat)

  depths = depth(f, fmat)
  f_obs = nrow(fmat)
  points = seq(1, f_obs) / f_obs

  cdf = sapply(points, function(x) sum(depths <= x) / f_obs)
  return(cdf)
}

#' @export
point_ED <- function(f1_cdf, f2_cdf) {
  # Calculate if function 1 or function 2 is more extreme
  # returns 1 if f1 is more extreme, -1 if f2 is more extreme, and 0 if equivalent
  rows = length(f1_cdf)

  for (x in 1:rows) {
    diff = f1_cdf[x] - f2_cdf[x]
    if (diff != 0.0) {
      return(sign(diff))
    }
  }

  return(0)
}

#' @export
ED <- function(fmat) {
  # Takes a matrix of functions (each column is a function) and returns there ED ordering
  # from deepest to shallowest

  # get dims
  obs = nrow(fmat)
  funcs = ncol(fmat)

  # find the depth CDF for each function
  dCDF = matrix(0, obs, funcs)
  for (f in 1:funcs) {
    dCDF[,f] = depth_CDF(fmat[,f], fmat)
  }

  # for each depth CDF (f1) count the number of depth CDFs (f2) that it's more extreme than
  edepth = rep(0, funcs)
  for (f1 in 1:funcs) {
    gt = 0
    for (f2 in 1:funcs) {
      gt = gt + (point_ED(dCDF[,f1], dCDF[,f2]) > 0)
    }
    edepth[f1] = gt
  }

  return(edepth / funcs)
}



