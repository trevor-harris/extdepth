#' @export
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

#' @export
edepth = function(g, fmat) {

  # save the dimensions for convenience
  obs = nrow(fmat)
  fns = ncol(fmat)

  # find the depths of g
  g_depth = depth(g, fmat)

  # find the depths of each f in fmat
  fmat_depth = matrix(0, obs, fns)
  for (j in 1:fns) {
    fmat_depth[,j] = depth(fmat[,j], fmat)
  }

  # get the allowed r values (for calculating the dCDF)
  r = sort(unique(c(g_depth, fmat_depth)))

  # find dCDF of g
  g_dcdf = sapply(r, function(x) sum(g_depth <= x))

  # find the dCDFs of each f in fmat
  fmat_dcdf = matrix(0, length(r), ncol(fmat_depth))
  for (j in 1:fns) {
    fmat_dcdf[,j] = sapply(r, function(x) sum(fmat_depth[,j] <= x))
  }

  # compare the dCDF of g against each dCDF of fmat
  ed = 0
  for (j in 1:fns) {
    for (i in 1:nrow(fmat_dcdf)) {
      if(g_dcdf[i] != fmat_dcdf[i,j]) {
        ed = ed + (fmat_dcdf[i,j] > g_dcdf[i])
        break
      }
    }
  }
  return(ed / fns)
}
