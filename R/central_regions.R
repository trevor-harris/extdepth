# Central regions

central_region <- function(extremal_depths, alpha = 0.05) {
  # Takes a list of extremal depths and an alpha level and returns the functions corresponding to
  # the lower and upper bounds. Also returns the median function

  # lower
  lower = which(edepths == alpha)

  # upper
  upper = which(edepths == 1-alpha)

  # median
  median = which.max(edepths)

  return(list(lower = lower,
              upper = upper,
              median = median))
}

extremal_boxplot <- function() {

}
