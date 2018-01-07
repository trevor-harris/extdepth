

Rcpp::cppFunction('int sgn(doublex) {
                    return (x < 0) ? -1 : (x > 0);
                  }')

Rcpp::cppFunction('int pwED(NumericVector f1_cdf, NumericVector f2_cdf) {
                  int rows = f1_cdf.length();
                  NumericVector fdiff(rows);

                  fdiff = sign(f1_cdf - f2_cdf);
                  for (int i = 0; i < rows; i++) {
                    if (fdiff(i) != 0) {
                      return fdiff(i);
                    };
                  };
                  
                  return 0;
                  }')


rm(list = ls())
gc()

library(microbenchmark)
library(extdepth)

n = 50
S = matrix(rnorm(n*1001), n, 1001)

f1 = S[,1]
f2 = S[,2]

cdf1 = dCDF(f1, S)
cdf2 = dCDF(f2, S)
# 
# dCDF(f1, S)
# depth_CDF(f1, S)

# cdf calc check
microbenchmark(dCDF(f1, S), times = 10)
microbenchmark(depth_CDF(f1, S), times = 10)

# pointwise ed check
t = 10000
microbenchmark(pointwise_ED(cdf1, cdf2), times = t, unit = "us")
microbenchmark(point_ED(cdf1, cdf2), times = t, unit = "us")
microbenchmark(pwED(cdf1, cdf2), times = t, unit = "us")

# full check
microbenchmark(ED_c(S), times = 10)
microbenchmark(ED(S), times = 10)


sum(ED_c(S) == ED(S))
