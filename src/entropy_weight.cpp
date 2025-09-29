#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector entropy_weight_cpp(NumericMatrix mat) {
  int n = mat.nrow();
  int m = mat.ncol();

  // Step 1: Normalise each column to [0,1]
  for (int j = 0; j < m; j++) {
    double min_val = mat(0, j);
    double max_val = mat(0, j);
    for (int i = 1; i < n; i++) {
      if (mat(i, j) < min_val) min_val = mat(i, j);
      if (mat(i, j) > max_val) max_val = mat(i, j);
    }
    double range = max_val - min_val;
    if (range > 0) {
      for (int i = 0; i < n; i++) {
        mat(i, j) = (mat(i, j) - min_val) / range;
      }
    } else {
      for (int i = 0; i < n; i++) {
        mat(i, j) = 0.0;
      }
    }
  }

  // Step 2: Compute proportions p_ij
  NumericMatrix P(n, m);
  for (int j = 0; j < m; j++) {
    double col_sum = 0.0;
    for (int i = 0; i < n; i++) col_sum += mat(i, j);
    if (col_sum > 0) {
      for (int i = 0; i < n; i++) P(i, j) = mat(i, j) / col_sum;
    }
  }

  // Step 3: Compute entropy for each column
  NumericVector E(m);
  double k = 1.0 / std::log(n);
  for (int j = 0; j < m; j++) {
    double sum_p_logp = 0.0;
    for (int i = 0; i < n; i++) {
      if (P(i, j) > 0) sum_p_logp += P(i, j) * std::log(P(i, j));
    }
    E[j] = -k * sum_p_logp;
  }

  // Step 4: Compute diversity degree
  NumericVector d(m);
  for (int j = 0; j < m; j++) d[j] = 1.0 - E[j];

  // Step 5: Compute weights
  double sum_d = 0.0;
  for (int j = 0; j < m; j++) sum_d += d[j];
  NumericVector w(m);
  for (int j = 0; j < m; j++) w[j] = d[j] / sum_d;

  return w;
}
