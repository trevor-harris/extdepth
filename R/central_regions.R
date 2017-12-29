# Central regions

central_region <- function(extremal_depths, alpha) {
  # Takes a list of extremal depths and an alpha level and returns the functions corresponding to
  # the lower and upper bounds

  # lower
  lower = which(edepths < alpha)
  lower = which.max(edepths[lower])

  # upper
  upper = which(edepths > 1-alpha)
  upper = which.min(edepths[upper])

  return(list(lower = lower, upper = upper))
}

