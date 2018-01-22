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
  NumericVector d1(obs);
  NumericVector d2(obs);

  for (int f1 = 0; f1 < fns; f1++) {
    d1 = depths(_, f1);
    edepths(f1) = 0;

    for (int f2 = 0; f2 < fns; f2++) {
      if (f1 == f2) {continue;};

      d2 = depths(_, f2);
      for (int x = 1; x <= obs; x++) {
        double r = x / double(obs);
        double cdf1 = sum(d1 <= r) / double(obs);
        double cdf2 = sum(d2 <= r) / double(obs);

        if (cdf1 != cdf2) {
          edepths(f1) += (cdf1 > cdf2);
          break;
        }
      }
    }
  }

  return edepths / fns;
};




// //' @export
// // [[Rcpp::export]]
// double dCDF_r(NumericVector depths, double r) {
//   return double(sum(depths <= r)) / depths.length();
// }
//
// //' @export
// // [[Rcpp::export]]
// NumericMatrix depths(NumericMatrix fmat) {
//   // get matrix dimension
//   int obs = fmat.nrow();
//   int fns = fmat.ncol();
//
//   NumericMatrix depth(obs, fns);
//   for (int f = 0; f < fns; f++) {
//     depth(_, f) = fast_depth(fmat(_, f), fmat);
//   }
//   return depth;
// }
