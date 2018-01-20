# fastr_ed = function(S) {
#   ncol.S = ncol(S)
#   nrow.S = nrow(S)
#
#   depths = matrix(0, nrow.S, ncol.S)
#   for (f in 1:ncol.S) {
#     depths[,f] = fast_depth(S[,f], S)
#   }
#
#   edepths = rep(0, ncol.S)
#   for (f1 in 1:ncol.S) {
#     gt = 0
#     for (f2 in 1:ncol.S) {
#       gt = gt + fast_compare(depths[,f1], depths[,f2])
#     }
#     edepths[f1] = gt / ncol.S
#   }
#   return(edepths)
# }
#
# # fed = fastr_ed(S)
#
# ncol.S = ncol(S)
# nrow.S = nrow(S)
#
# depths = matrix(0, nrow.S, ncol.S)
# for (f in 1:ncol.S) {
#   depths[,f] = fast_depth(S[,f], S)
# }
#
# cdf1 = dCDF(depths[,1])
# cdf2 = sapply(1:nrow.S, function(x) dCDF_r(depths[,1], x/nrow.S))


fastr_ed = function(S) {
  ncol.S = ncol(S)
  nrow.S = nrow(S)

  depths = matrix(0, nrow.S, ncol.S)
  for (f in 1:ncol.S) {
    depths[,f] = fast_depth(S[,f], S)
  }

  edepths = rep(0, ncol.S)
  for (f1 in 1:ncol.S) {
    gt = 0
    d1 = depths[,f1]
    for (f2 in 1:ncol.S) {
      if (f1 == f2) next

      d2 = depths[,f2]
      for (x in 1:nrow.S) {
        r = x/nrow.S
        cdf1 = dCDF_r(d1, r)
        cdf2 = dCDF_r(d2, r)
        if (cdf1 != cdf2) {
          gt = gt + (cdf1 > cdf2)
          break
        }
      }
    }
    edepths[f1] = gt / ncol.S
  }
  return(edepths)
}
fastr_ed = compiler::cmpfun(fastr_ed)

microbenchmark(ED(S), times = 10, unit = "s")
microbenchmark(fastr_ed(S), times = 10, unit = "s")
microbenchmark(EDr(S), times = 10, unit = "s")







