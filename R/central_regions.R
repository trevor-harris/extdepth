# Central regions

#' @export
central_region <- function(fmat, ext.depths, alpha = 0.05) {
  # Takes a list of extremal depths and an alpha level and returns the functions corresponding to
  # the lower and upper bounds. Also returns the median function

  # filter out functions outside the alpha level
  fset = fmat[,ext.depths >= alpha]

  # lower
  lower = sapply(1:nrow(fset), function(x) min(fset[x,]))

  # upper
  upper = sapply(1:nrow(fset), function(x) max(fset[x,]))

  return(list(lower = lower,
              upper = upper))
}

extremal_boxplot <- function() {

}
