# Central regions

central_region <- function(ext.depths, alpha = 0.05) {
  # Takes a list of extremal depths and an alpha level and returns the functions corresponding to
  # the lower and upper bounds. Also returns the median function

  # lower
  lower = which(ext.depths == alpha)

  # upper
  upper = which(ext.depths == 1-alpha)

  # median
  median = which.max(ext.depths)

  return(list(lower = lower,
              upper = upper,
              median = median))
}

extremal_boxplot <- function() {

}
