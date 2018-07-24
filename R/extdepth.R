
#' Compute the extremal depth of a function \code{g} with respect to a set of functions \code{fmat}
#'
#' @param g Vector. Function that you want to find the extremal depth of.
#' @param fmat Matrix of functions. Each column is a function.
#' @param depth_function either "standard" (default) or "rank". Standard orders functions from the center
#'  outward and rank orders least to greatest.
#'
#' @return The extremal depth value as a real number.
#'
#' @note \code{g} and \code{fmat} should have the same number of observations for each function. Furthermore
#' functions should be sampled at the same time points.
#' @export
edepth = function(g, fmat, depth_function = "standard") {

  # save the dimensions for convenience
  obs = nrow(fmat)
  fns = ncol(fmat)

  # find the depths of g
  g_depth = depth(g, fmat)

  # find the depths of each f in fmat
  if (depth_function == "standard") fmat_depth = depth_set(fmat)
  if (depth_function == "rank") fmat_depth = rank_depth(fmat)

  # get the allowed r values (for calculating the dCDF)
  rvals = sort(unique(c(g_depth, fmat_depth)))

  partial_ed = rep(1, fns)
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

#' Compute the extremal depths for all functions in \code{fmat}, with respect to \code{fmat}
#'
#' @param fmat Matrix of functions. Each column is a function.
#' @param depth_function either "standard" (default) or "rank". Standard orders functions from the center
#'  outward and rank orders least to greatest.
#'
#' @return A vector of extremal depth values that correspond to the functions in fmat.
#' @export
edepth_set = function(fmat, depth_function = "standard") {
  # find the depths of each f in fmat
  if (depth_function == "standard") depths = depth_set(fmat)
  if (depth_function == "rank") depths = rank_depth(fmat)

  arr = 1:ncol(depths)
  rvals = unique(sort(depths))
  ed = order(func_quickSort(depths, arr, rvals)) / ncol(fmat)

  return(1-ed)
}

#' Compute the extremal depths for all functions in \code{fmat}, with respect to \code{fmat}, then orders them from least to
#' greatest extremal depth value.
#'
#' @param fmat Matrix of functions. Each column is a function.
#' @param depth_function either "standard" (default) or "rank". Standard orders functions from the center
#'  outward and rank orders least to greatest.
#'
#' @return A vector of extremal depth values that correspond to the functions in fmat.
#' @export
edepth_sort = function(fmat, depth_function = "standard") {
  ed.set = edepth_set(fmat, depth_function)
  ranks = match(sort(ed.set, decreasing = F), ed.set)
  return(fmat[,ranks])
}


# edepth_set = function(fmat, depth_function = "standard") {
#
#   # save the dimensions for convenience
#   obs = nrow(fmat)
#   fns = ncol(fmat)
#
#   # find the depths of each f in fmat
#   if (depth_function == "standard") fmat_depth = depth_set(fmat)
#   if (depth_function == "rank") fmat_depth = rank_depth(fmat)
#
#   # get the allowed r values (for calculating the dCDF)
#   rvals = unique(sort(fmat_depth))
#
#   # only need to fill in the top triangle
#   # invert later to get full partial_eds
#   partial_ed = matrix(1, fns, fns)
#
#   for (f1 in 1:(fns-1)) {
#     for (f2 in (f1+1):fns) {
#       for (r in rvals) {
#         d1 = sum(fmat_depth[,f1] <= r)
#         d2 = sum(fmat_depth[,f2] <= r)
#         if(d1 != d2) {
#           partial_ed[f1, f2] = (d2 > d1)
#           partial_ed[f2, f1] = (d2 < d1)
#           break;
#         }
#       }
#     }
#   }
#
#   ed = rowMeans(partial_ed)
#   return(ed)
# }

# edepth_multi = function(gmat, sorted_fmat, fdepths) {
#
#   # shorten name
#   fmat = sorted_fmat
#
#   # save the dimensions for convenience
#   g.obs = nrow(gmat)
#   g.fns = ncol(gmat)
#
#   f.obs = nrow(fmat)
#   f.fns = ncol(fmat)
#
#   # find the depths of each f in fmat
#   fmat_depth = fdepths
#
#   # find the depths of each g in gmat
#   gmat_depth = depth_set(gmat)
#
#   # get the allowed r values (for calculating the dCDF)
#   rvals = unique(sort(c(gmat_depth, fmat_depth)))
#
#   partial_ed = matrix(0, g.fns, f.fns)
#   for (g in 1:g.fns) {
#     for (f in 1:f.fns) {
#       for (r in rvals) {
#         d1 = sum(gmat_depth[,g] <= r)
#         d2 = sum(fmat_depth[,f] <= r)
#         if(d1 != d2) {
#           partial_ed[g, f] = (d2 > d1)
#           break;
#         }
#       }
#       if (partial_ed[g, f]==0) break;
#     }
#   }
#   return(rowSums(partial_ed) / f.fns)
# }
