# Central regions

#' Find the upper and lower bound for the \code{1-alpha} central region for a set of functions \code{fmat}
#'
#' @param fmat Matrix of functions. Each column is a function.
#' @param edepths Vector of extremal depths for the functions in fmat. first value should correspond
#' to the first function etc.
#' @param alpha Numeric. Default is 0.05.
#'
#' @return A list containing the upper and lower bound as separate vectors.
#' @export
central_region <- function(fmat, edepths, alpha = 0.05) {
  # Takes a list of extremal depths and an alpha level and returns the functions corresponding to
  # the lower and upper bounds. Also returns the median function

  # filter out functions outside the alpha level
  fset = fmat[,edepths > alpha]

  # lower
  lower = sapply(1:nrow(fset), function(x) min(fset[x,]))

  # upper
  upper = sapply(1:nrow(fset), function(x) max(fset[x,]))

  return(list(lower = lower,
              upper = upper))
}

extremal_boxplot <- function() {

}
