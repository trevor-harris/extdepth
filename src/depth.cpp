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

  for (int row = 0; row < rows; row++) {
    depths(row) = 1 - (double(std::abs(sum(sign(f(row) - fmat(row, _))))) / double(cols));
  };
  return depths;
}

//' @export
// [[Rcpp::export]]
NumericVector extdepth(NumericVector g, NumericMatrix fmat) {

  // get matrix dimension
  int obs = fmat.nrow();
  int fns = fmat.ncol();

  // compute the depth of G wrt to fmat
  NumericVector g_depth = depth(g, fmat);

  // compute the depths for each function in fmat
  NumericMatrix depths(obs, fns);
  for (int f = 0; f < fns; f++) {
    depths(_, f) = depth(fmat(_, f), fmat);
  };

  // compute number of functions each function is greater than
  double edepth = 0;
  double cdf1 = 0;
  double cdf2 = 0;
  double r = 0;

  for (int f = 0; f < fns; f++) {

    NumericMatrix::Column f_depth = depths(_, f);
    for (int x = 1; x <= fns; x++) {
      r = x / double(fns);
      cdf1 = sum(g_depth <= r);
      cdf2 = sum(f_depth <= r);

      if (cdf1 != cdf2) {
        edepth += (cdf1 > cdf2);
        break;
      }
    }
  }

  return edepth / fns;
};

//' @export
// [[Rcpp::export]]
NumericVector extdepth_all(NumericMatrix fmat) {

  // get matrix dimension
  int obs = fmat.nrow();
  int fns = fmat.ncol();

  // compute the dCDF for each function in fmat
  NumericMatrix depths(obs, fns);
  for (int f = 0; f < fns; f++) {
    depths(_, f) = depth(fmat(_, f), fmat);
  };

  // compute number of functions each function is greater than
  NumericVector edepths(fns);
  double cdf1 = 0;
  double cdf2 = 0;
  double r = 0;

  for (int f1 = 0; f1 < fns; f1++) {
    edepths(f1) = 1;

    NumericMatrix::Column d1 = depths(_, f1);
    for (int f2 = 0; f2 < fns; f2++) {
      if (f1 == f2) {continue;};

      NumericMatrix::Column d2 = depths(_, f2);
      for (int x = 1; x <= fns; x++) {
        r = x / double(fns);
        cdf1 = sum(d1 <= r);
        cdf2 = sum(d2 <= r);

        if (cdf1 != cdf2) {
          edepths(f1) += (cdf1 > cdf2);
          break;
        }
      }
    }
  }

  return edepths / fns;
};

