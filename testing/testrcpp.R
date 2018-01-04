# Rcpp::cppFunction('double dCDF(NumericVector f, NumericMatrix fmat) {
# 
#                   for (int row = 0; row < rown; row++) {
#                   double r = (row + 1) / rown;
#                   };
#                   
#                   }')

Rcpp::cppFunction('NumericVector dCDF1(NumericVector f, NumericMatrix fmat) {

                  // get matrix dimension
                  int rown = fmat.nrow();
                  int coln = fmat.ncol();
                  
                  // compute the depths at each observation of f
                  NumericVector depths(rown);
                  LogicalVector diffs(coln);

                  for(int row = 0; row < rown; row++) {
                    diffs = (f(row) > fmat(row, _)) - (f(row) < fmat(row, _));
                    depths(row) = 1 - (std::abs(sum(diffs)) / coln);
                  };
                  
                    // compute the cdf
                  NumericVector cdf(rown);
                  LogicalVector lt_r(rown);
                  double r;
                  for (int r_ind = 0; r_ind < rown; r_ind++) {
                  r = double((r_ind + 1)) / double(rown);
                  
                  lt_r = depths <= r;
                  cdf(r_ind) = double(sum(lt_r)) / double(rown);
                  };
                  
                  return cdf;
                  }')



rm(list = ls())
gc()

library(microbenchmark)
library(extdepth)

n = 5000
S = matrix(rnorm(n*n), n, 1000)

f1 = S[,1]
f2 = S[,2]

cdf1 = dCDF(f1, S)
cdf2 = dCDF(f2, S)

dCDF(f1, S)
depth_CDF(f1, S)

# cdf calc check
microbenchmark(dCDF(f1, S), times = 100)
microbenchmark(depth_CDF(f1, S), times = 100)

# pointwise ed check
microbenchmark(pointwise_ED(cdf1, cdf2), times = 200)
microbenchmark(point_ED(cdf1, cdf2), times = 200)

# full check
microbenchmark(ED_c(S), times = 1)
microbenchmark(ED(S), times = 1)




