#include <Rcpp.h>
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
NumericVector fast_depth(NumericVector f, NumericMatrix fmat) {

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
double dCDF_r(NumericVector depths, double r) {
  return double(sum(depths <= r)) / depths.length();
}

// make this do the very innermost loop
//' @export
// [[Rcpp::export]]
int fast_compare(NumericVector d1, NumericVector d2) {
  int rows = d1.length();
  double r;
  double cdf_diff;

  // left tail comparions of the cdfs. returns 1 if f1 > f2 and -1 if f1 < f2.
  for (int row = 1; row < rows+1; row++) {
    r = row / double(rows);
    cdf_diff = (sum(d1 <= r) - sum(d2 <= r)) / double(rows);
    //cdf_diff = dCDF_r(d1, r) - dCDF_r(d2, r);
    if (cdf_diff > 0) {
      return 1;
    };
  };

  // if cdfs are identical then return 0. f1 ~ f2.
  return 0;
}



//' @export
// [[Rcpp::export]]
NumericVector fast_ED(NumericMatrix fmat) {

  // get matrix dimension
  int obs = fmat.nrow();
  int fns = fmat.ncol();

  // compute the dCDF for each function in fmat
  NumericMatrix depths(obs, fns);
  for (int f = 0; f < fns; f++) {
    depths(_, f) = fast_depth(fmat(_, f), fmat);
  };

  // compute number of functions each function is greater than
  NumericVector edepths(fns);
  for (int f1 = 0; f1 < fns; f1++) {
    for (int f2 = 0; f2 < fns; f2++) {
      edepths(f1) += fast_compare(depths(_, f1), depths(_, f2));
    };
  };

  return edepths / fns;
};

