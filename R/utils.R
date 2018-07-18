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

ed_compare = function(f, m, depths, rvals) {
  for (r in rvals) {
    df = sum(depths[,f] <= r)
    dm = sum(depths[,m] <= r)
    if(df != dm) {
      if((df > dm)) return(FALSE)
      else return(TRUE)
    }
  }
}

# quickly sort depths
func_quickSort <- function(depths, arr, rvals) {
  # Pick a number at random.
  mid <- sample(arr, 1)

  # Place-holders for left and right values.
  left <- c()
  right <- c()

  # Move all the smaller values to the left, bigger values to the right.
  lapply(arr[arr != mid], function(d) {
    if (ed_compare(d, mid, depths, rvals)) {
      left <<- c(left, d)
    }
    else {
      right <<- c(right, d)
    }
  })

  if (length(left) > 1) {
    left <- func_quickSort(depths, left, rvals)
  }

  if (length(right) > 1) {
    right <- func_quickSort(depths, right, rvals)
  }

  # Finally, return the sorted values.
  c(left, mid, right)
}

