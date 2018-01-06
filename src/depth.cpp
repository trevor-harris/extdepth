#include <Rcpp.h>
using namespace Rcpp;


// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//' @export
// [[Rcpp::export]]
NumericVector depth_c(NumericVector f, NumericMatrix fmat) {

  // get matrix dimension
  int rows = fmat.nrow();
  int cols = fmat.ncol();

  // compute the depths at each observation of f
  NumericVector depth(rows);

  for(int row = 0; row < rows; row++) {
    depth(row) = 1 - (double(std::abs(sum(sign(f(row) - fmat(row, _))))) / double(cols));
  };
  return depth;
}

//' @export
// [[Rcpp::export]]
NumericVector dCDF(NumericVector f, NumericMatrix fmat) {

  // get matrix dimension
  int rows = fmat.nrow();
  int cols = fmat.ncol();

  // compute the depths at each observation of f
  NumericVector depths(rows);

  for(int row = 0; row < rows; row++) {
    depths(row) = 1 - (double(std::abs(sum(sign(f(row) - fmat(row, _))))) / double(cols));
  };

  // compute the cdf
  NumericVector cdf(rows);
  double r;

  for (int r_ind = 0; r_ind < rows; r_ind++) {
    r = double(r_ind + 1) / double(rows);
    cdf(r_ind) = double(sum(depths <= r)) / double(rows);
  };

  return cdf;
}

//' @export
// [[Rcpp::export]]
int pointwise_ED(NumericVector f1_cdf, NumericVector f2_cdf) {
  int rows = f1_cdf.length();
  double differ;

  for (int row = 0; row < rows; row++) {
    differ = f1_cdf(row) - f2_cdf(row);
    if (differ != 0.0) {
      return (differ > 0) - (differ < 0);
    };
  };

  return 0;
}

//' @export
// [[Rcpp::export]]
NumericVector ED_c(NumericMatrix fmat) {

  // get matrix dimension
  int obs = fmat.nrow();
  int fns = fmat.ncol();

  // compute the dCDF for each function in fmat
  NumericMatrix cdf(obs, fns);
  for (int f = 0; f < fns; f++) {
    cdf(_, f) = dCDF(fmat(_, f), fmat);
  };

  // compute number of functions each function is greater than
  NumericVector edepths(fns);
  for (int f1 = 0; f1 < fns; f1++) {
    for (int f2 = 0; f2 < fns; f2++) {
      edepths(f1) += (pointwise_ED(cdf(_, f1), cdf(_, f2)) > 0);
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
