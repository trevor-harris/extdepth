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
NumericVector dCDF(NumericVector f, NumericMatrix fmat) {

  // get matrix dimension
  int rows = fmat.nrow();
  int cols = fmat.ncol();

  // compute the depths at each observation of f
  NumericVector depths(rows);
  // NumericVector diffs(cols);

  for(int row = 0; row < rows; row++) {
    // diffs = (f(row) > fmat(row, _)) - (f(row) < fmat(row, _));
    // diffs = sign(f(row) - fmat(row, _));
    depths(row) = 1 - (double(std::abs(sum(sign(f(row) - fmat(row, _))))) / double(cols));
  };

  // compute the cdf
  NumericVector cdf(rows);
  // LogicalVector lt_r(rows);
  double r;

  for (int r_ind = 0; r_ind < rows; r_ind++) {
    r = double(r_ind + 1) / double(rows);

    // lt_r = depths <= r;
    cdf(r_ind) = double(sum(depths <= r)) / double(rows);
  };

  return cdf;
}

//' @export
// [[Rcpp::export]]
int pointwise_ED(NumericVector f1_cdf, NumericVector f2_cdf) {
  int rows = f1_cdf.length();
  double cdf_diff;

  for (int row = 0; row < rows; row++) {
    cdf_diff = f1_cdf(row) - f2_cdf(row);
    if (cdf_diff != 0.0) {
      return (cdf_diff > 0) - (cdf_diff < 0);
    };
  };

  return 0;
}

//' @export
// [[Rcpp::export]]
NumericVector ED_c(NumericMatrix fmat) {

  // get matrix dimension
  int rows = fmat.nrow();
  int cols = fmat.ncol();

  NumericMatrix cdf(rows, cols);
  for (int col = 0; col < cols; col++) {
    cdf(_, col) = dCDF(fmat(_, col), fmat);
  };

  NumericVector edepths(cols);
  int gt;
  for (int col1 = 0; col1 < cols; col1++) {
    gt = 0;
    for (int col2 = 0; col2 < cols; col2++) {
        gt += (pointwise_ED(cdf(_, col1), cdf(_, col2)) > 0);
    };

    edepths(col1) = gt;
  };

  return edepths / rows;
};





// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
*/
