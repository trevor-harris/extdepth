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
  rvals = sort(unique(c(g_depth, fmat_depth)))

  partial_ed = rep(0, fns)
  for (f in 1:fns) {
    for (r in rvals) {
      dg = sum(g_depth <= r)
      df = sum(fmat_depth[,f] <= r)
      if(dg != df) {
        partial_ed[f] = (df > dg)
        break;
      }
    }
  }
  ed = mean(partial_ed)
  return(ed)
}


#' @export
edepth_multi_fast = function(gmat, sorted_fmat, fdepths) {

  # shorten name
  fmat = sorted_fmat

  # save the dimensions for convenience
  g.obs = nrow(gmat)
  g.fns = ncol(gmat)

  f.obs = nrow(fmat)
  f.fns = ncol(fmat)

  # find the depths of each f in fmat
  fmat_depth = fdepths

  # find the depths of each f in fmat
  gmat_depth = matrix(0, g.obs, g.fns)
  for (j in 1:g.fns) {
    gmat_depth[,j] = depth(gmat[,j], fmat)
  }

  # get the allowed r values (for calculating the dCDF)
  rvals = unique(sort(c(gmat_depth, fmat_depth)))

  partial_ed = matrix(0, g.fns, f.fns)
  for (g in 1:g.fns) {
    for (f in 1:f.fns) {
      for (r in rvals) {
        d1 = sum(gmat_depth[,g] <= r)
        d2 = sum(fmat_depth[,f] <= r)
        if(d1 != d2) {
          partial_ed[g, f] = (d2 > d1)
          break;
        }
      }
      if (partial_ed[g, f]==0) break;
    }
  }
  return(rowSums(partial_ed) / f.fns)
}


#' @export
edepth_set = function(fmat) {

  # save the dimensions for convenience
  obs = nrow(fmat)
  fns = ncol(fmat)

  # find the depths of each f in fmat
  fmat_depth = matrix(0, obs, fns)
  for (j in 1:fns) {
    fmat_depth[,j] = depth(fmat[,j], fmat)
  }

  # get the allowed r values (for calculating the dCDF)
  rvals = unique(sort(fmat_depth))

  # only need to fill in the top triangle
  # invert later to get full partial_eds
  partial_ed_upper = matrix(0, fns, fns)
  for (f1 in 1:(fns-1)) {
    for (f2 in (f1+1):fns) {
      for (r in rvals) {
        d1 = sum(fmat_depth[,f1] <= r)
        d2 = sum(fmat_depth[,f2] <= r)
        if(d1 != d2) {
          partial_ed_upper[f1, f2] = (d2 > d1)
          break;
        }
      }
    }
  }

  # tranpose and invert the top triangle
  partial_ed_lower = abs(1 - t(partial_ed_upper))*lower.tri(partial_ed_upper)
  partial_ed = partial_ed_lower + partial_ed_upper + diag(1, fns, fns)
  ed = rowMeans(partial_ed)

  return(ed)
}

