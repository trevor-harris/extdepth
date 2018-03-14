#' @export
sort_by_ed = function(fmat) {
  ed.set = edepth_set(fmat)
  ranks = match(sort(ed.set, decreasing = F), ed.set)
  return(fmat[,ranks])
}