// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <algorithm>
#include <cmath>

using namespace Rcpp;
using namespace RcppParallel;

// ----- Parallel Worker for TSI -----
struct TSIWorker : public Worker {
  const RMatrix<double> mat;
  RVector<double> tsi;

  TSIWorker(const NumericMatrix mat, NumericVector tsi)
    : mat(mat), tsi(tsi) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      double maxVal = R_NegInf;
      double sumVal = 0.0;
      int n = mat.ncol();

      for (int j = 0; j < n; j++) {
        double v = mat(i, j);
        if (v > maxVal) maxVal = v;
      }
      if (!R_finite(maxVal) || maxVal <= 0.0) {
        tsi[i] = NA_REAL;
        continue;
      }
      for (int j = 0; j < n; j++) {
        sumVal += 1.0 - (mat(i, j) / maxVal);
      }
      tsi[i] = sumVal / (n - 1);
    }
  }
};

// ----- Parallel Worker for TSE -----
struct TSEWorker : public Worker {
  const RMatrix<double> mat;
  RVector<double> tse;

  TSEWorker(const NumericMatrix mat, NumericVector tse)
    : mat(mat), tse(tse) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      double maxVal = R_NegInf;
      double minVal = R_PosInf;
      double sumVal = 0.0;
      int n = mat.ncol();

      for (int j = 0; j < n; j++) {
        double v = mat(i, j);
        if (v > maxVal) maxVal = v;
        if (v < minVal) minVal = v;
      }
      if (!R_finite(maxVal) || maxVal <= 0.0) {
        tse[i] = NA_REAL;
        continue;
      }
      for (int j = 0; j < n; j++) {
        sumVal += 1.0 - (minVal / maxVal);
      }
      tse[i] = sumVal / (n - 1);
    }
  }
};

// [[Rcpp::export]]
List compute_TSI_TSE_parallel(NumericMatrix expr) {
  int nGenes = expr.nrow();
  NumericVector tsi(nGenes), tse(nGenes);

  TSIWorker tsi_worker(expr, tsi);
  TSEWorker tse_worker(expr, tse);

  parallelFor(0, nGenes, tsi_worker);
  parallelFor(0, nGenes, tse_worker);

  return List::create(
    _["TSI"] = tsi,
    _["TSE"] = tse
  );
}
