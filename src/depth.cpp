#include <Rcpp.h>
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
NumericVector depth(NumericVector f, NumericMatrix fmat) {

  // get matrix dimension
  int rows = fmat.nrow();
  int cols = fmat.ncol();

  // compute the depths at each observation of f
  NumericVector depths(rows);

  for(int row = 0; row < rows; row++) {
    depths(row) = 1 - (double(std::abs(sum(sign(f(row) - fmat(row, _))))) / double(cols));
  };
  return depths;
}

//' @export
// [[Rcpp::export]]
NumericVector dCDF(NumericVector depths) {

  // get matrix dimension
  double obs = depths.length();

  // compute the cdf
  NumericVector cdf(obs);
  double r;

  for (int r_ind = 0; r_ind < obs; r_ind++) {
    r = double(r_ind + 1) / double(obs);
    cdf(r_ind) = double(sum(depths <= r)) / double(obs);
  };

  return cdf;
}

//' @export
// [[Rcpp::export]]
int ed_compare(NumericVector f1_cdf, NumericVector f2_cdf) {
  int rows = f1_cdf.length();
  double cdf_diff;

  // left tail comparions of the cdfs. returns 1 if f1 > f2 and -1 if f1 < f2.
  for (int row = 0; row < rows; row++) {
    cdf_diff = f1_cdf(row) - f2_cdf(row);
    if (cdf_diff != 0.0) {
      return (cdf_diff > 0) - (cdf_diff < 0);
    };
  };

  // if cdfs are identical then return 0. f1 ~ f2.
  return 0;
}

//' @export
// [[Rcpp::export]]
NumericVector ED(NumericMatrix fmat) {

  // get matrix dimension
  int obs = fmat.nrow();
  int fns = fmat.ncol();

  // compute the dCDF for each function in fmat
  NumericMatrix cdf(obs, fns);
  for (int f = 0; f < fns; f++) {
    cdf(_, f) = dCDF(depth(fmat(_, f), fmat));
  };

  // compute number of functions each function is greater than
  NumericVector edepths(fns);
  for (int f1 = 0; f1 < fns; f1++) {
    for (int f2 = 0; f2 < fns; f2++) {
      edepths(f1) += (ed_compare(cdf(_, f1), cdf(_, f2)) > 0);
    };
  };

  return edepths / fns;
};





// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
*/
