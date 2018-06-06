# compute the pointwise depth, with respect to fmat, for the function g.
depth <- function(g, fmat) {

  # Computes the depth values of a function with respect to a set of functions (fmat)
  fn = ncol(fmat)
  depth = rep(0, length(g))

  for (row in 1:nrow(fmat)) {
    diff = abs(sum(sign(g[row] - fmat[row,])))
    depth[row] = 1 - (diff / fn)
  }

  return(depth)
}

# compute the pointwise depth, with respect to fmat, for each function in fmat
#' @export
depth_set <- function(fmat) {
  # save the dimensions for convenience
  obs = nrow(fmat)
  fns = ncol(fmat)

  # find the depths of each f in fmat
  fmat_depth = matrix(0, obs, fns)
  for (j in 1:fns) {
    fmat_depth[,j] = depth(fmat[,j], fmat)
  }
  return(fmat_depth)
}

# compute the pointwise ranking, with respect to fmat, for each function in fmat
rank_depth = function(fmat) {
  # technically not a depth
  # At each time t, order the function observations from least to greatest
  obs = nrow(fmat)
  fns = ncol(fmat)

  rdepth = matrix(0, obs, fns)
  for (i in 1:obs) {
    rdepth[i,] = (rank(fmat[i,])-1) / fns
  }
  return(1 - rdepth)
}
