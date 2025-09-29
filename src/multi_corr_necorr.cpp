// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <algorithm>
#include <cmath>
#include <random>
#include <vector>
#include <numeric>
#include <utility>

using namespace Rcpp;
using namespace RcppParallel;
using namespace std;

// ===================== Pearson correlation =====================
inline double pearsonCorrelation(const vector<double> &xData, const vector<double> &yData) {
  int n = xData.size();
  double mean_x = std::accumulate(xData.begin(), xData.end(), 0.0) / n;
  double mean_y = std::accumulate(yData.begin(), yData.end(), 0.0) / n;
  double num = 0.0, den_x = 0.0, den_y = 0.0;
  for (int i = 0; i < n; i++) {
    double dx = xData[i] - mean_x;
    double dy = yData[i] - mean_y;
    num += dx * dy;
    den_x += dx * dx;
    den_y += dy * dy;
  }
  if (den_x <= 0.0 || den_y <= 0.0) return 0.0;
  return num / sqrt(den_x * den_y);
}

// ===================== Spearman correlation =====================
inline double spearmanCorrelation(const vector<int> &xRank, const vector<int> &yRank) {
  int n = xRank.size();
  double mean_x = (n - 1) / 2.0;
  double mean_y = (n - 1) / 2.0;
  double num = 0.0, den_x = 0.0, den_y = 0.0;
  for (int i = 0; i < n; i++) {
    double dx = xRank[i] - mean_x;
    double dy = yRank[i] - mean_y;
    num += dx * dy;
    den_x += dx * dx;
    den_y += dy * dy;
  }
  if (den_x <= 0.0 || den_y <= 0.0) return 0.0;
  return num / sqrt(den_x * den_y);
}

// ===================== Kendall correlation =====================
inline double kendallCorrelation(const vector<double> &xData, const vector<double> &yData) {
  int n = xData.size();
  int concordant = 0, discordant = 0;
  for (int i = 0; i < n - 1; i++) {
    for (int j = i + 1; j < n; j++) {
      double dx = xData[i] - xData[j];
      double dy = yData[i] - yData[j];
      double prod = dx * dy;
      if (prod > 0) concordant++;
      else if (prod < 0) discordant++;
    }
  }
  double denom = 0.5 * n * (n - 1);
  if (denom == 0) return 0.0;
  return (concordant - discordant) / denom;
}

// ===================== Empirical CDF helper =====================
vector<double> empiricalCDF(const vector<double> &data) {
  int n = data.size();
  vector<pair<double, int>> sorted(n);
  for (int i = 0; i < n; i++) {
    sorted[i] = {data[i], i};
  }
  sort(sorted.begin(), sorted.end());
  vector<double> cdf(n);
  for (int rank = 0; rank < n; rank++) {
    cdf[sorted[rank].second] = (rank + 1.0) / n;
  }
  return cdf;
}

// ===================== Asymmetric Gini correlation γ(X,Y) =====================
double giniAsym(const vector<double> &X, const vector<double> &Y) {
  int n = X.size();
  double minX = *min_element(X.begin(), X.end());
  double maxX = *max_element(X.begin(), X.end());
  if (minX == maxX) return 0.0;
  double minY = *min_element(Y.begin(), Y.end());
  double maxY = *max_element(Y.begin(), Y.end());
  if (minY == maxY) return 0.0;

  vector<double> Fy = empiricalCDF(Y);
  vector<double> Fx = empiricalCDF(X);

  double meanX = accumulate(X.begin(), X.end(), 0.0) / n;
  double meanFy = accumulate(Fy.begin(), Fy.end(), 0.0) / n;
  double cov_X_Fy = 0.0;
  for (int i = 0; i < n; i++) {
    cov_X_Fy += (X[i] - meanX) * (Fy[i] - meanFy);
  }

  double meanFx = accumulate(Fx.begin(), Fx.end(), 0.0) / n;
  double cov_X_Fx = 0.0;
  for (int i = 0; i < n; i++) {
    cov_X_Fx += (X[i] - meanX) * (Fx[i] - meanFx);
  }

  if (cov_X_Fx == 0.0) return 0.0;
  return cov_X_Fy / cov_X_Fx;
}

// ===================== Symmetric Gini correlation (Ma et al., 2012) =====================
double calcGiniSymmetric(const vector<double> &X, const vector<double> &Y) {
  double gammaXY = giniAsym(X, Y);
  double gammaYX = giniAsym(Y, X);
  if (!R_finite(gammaXY) || !R_finite(gammaYX)) return 0.0;

  vector<double> rankY = empiricalCDF(Y);
  int n = X.size();
  double meanX = accumulate(X.begin(), X.end(), 0.0) / n;
  double meanRankY = accumulate(rankY.begin(), rankY.end(), 0.0) / n;
  double cov_X_rankY = 0.0;
  for (int i = 0; i < n; i++) {
    cov_X_rankY += (X[i] - meanX) * (rankY[i] - meanRankY);
  }

  double signVal = (cov_X_rankY > 0) ? 1.0 : ((cov_X_rankY < 0) ? -1.0 : 0.0);
  return signVal * sqrt(fabs(gammaXY * gammaYX));
}

// ===================== Optimized permutation test =====================
double permutationTests(double originalGCC, const vector<double> &xData, const vector<int> &xIdx,
                        const vector<double> &yData, const vector<int> &yIdx,
                        const vector<double> &wt, int perm) {
  int n = xData.size();
  if (fabs(originalGCC) >= 0.999999) {
    return 1.0 / (perm + 1); // extremely significant
  }

  double xDenominator = 0, yDenominator = 0;
  for (int i = 0; i < n; i++) {
    xDenominator += wt[i] * xData[xIdx[i]];
    yDenominator += wt[i] * yData[yIdx[i]];
  }
  if (xDenominator == 0 || yDenominator == 0) return 1.0;

  vector<int> indices(n);
  iota(indices.begin(), indices.end(), 0);
  mt19937 rng(random_device{}());

  int m = 0;
  for (int p = 0; p < perm; p++) {
    shuffle(indices.begin(), indices.end(), rng);

    double xNumerator = 0, yNumerator = 0;
    for (int i = 0; i < n; i++) {
      xNumerator += wt[i] * xData[indices[yIdx[i]]];
      yNumerator += wt[i] * yData[indices[xIdx[i]]];
    }

    double xGCC = xNumerator / xDenominator;
    double yGCC = yNumerator / yDenominator;

    if (fabs(xGCC) >= fabs(originalGCC) || fabs(yGCC) >= fabs(originalGCC)) {
      m++;
    }

    if ((p == 10 && ((double)m / p) > 0.5) ||
        (p == 100 && ((double)m / p) > 0.2) ||
        (p == 500 && ((double)m / p) > 0.1)) {
      perm = p + 1;
      break;
    }
  }

  return ((double)m + 1) / ((double)perm + 1); // +1 for stability
}

// ===================== Parallel Worker =====================
struct MultiCorrWorker : public Worker {
  const RMatrix<double> data;
  const RMatrix<int> ranks;
  const RVector<double> wt;
  const RVector<int> src;
  const RVector<int> tgt;
  int nSamples;
  int bootstrapIterations;
  bool useBestGCC;
  bool asymmetricGCC;

  RVector<double> gcc, gcc_pval, gcc1_vec, gcc2_vec;
  RVector<double> pcc, pcc_pval, scc, scc_pval, kcc, kcc_pval;

  MultiCorrWorker(const NumericMatrix data,
                  const IntegerMatrix ranks,
                  const NumericVector wt,
                  const IntegerVector src,
                  const IntegerVector tgt,
                  int nSamples,
                  int bootstrapIterations,
                  bool useBestGCC,
                  bool asymmetricGCC,
                  NumericVector gcc, NumericVector gcc_pval,
                  NumericVector gcc1_vec, NumericVector gcc2_vec,
                  NumericVector pcc, NumericVector pcc_pval,
                  NumericVector scc, NumericVector scc_pval,
                  NumericVector kcc, NumericVector kcc_pval)
    : data(data), ranks(ranks), wt(wt),
      src(src), tgt(tgt),
      nSamples(nSamples),
      bootstrapIterations(bootstrapIterations),
      useBestGCC(useBestGCC),
      asymmetricGCC(asymmetricGCC),
      gcc(gcc), gcc_pval(gcc_pval),
      gcc1_vec(gcc1_vec), gcc2_vec(gcc2_vec),
      pcc(pcc), pcc_pval(pcc_pval),
      scc(scc), scc_pval(scc_pval),
      kcc(kcc), kcc_pval(kcc_pval) {}

  void operator()(size_t begin, size_t end) {
    vector<double> xData(nSamples), yData(nSamples), wtVec(nSamples);
    vector<int> xRank(nSamples), yRank(nSamples);
    wtVec.assign(wt.begin(), wt.end());

    for (size_t i = begin; i < end; i++) {
      int i1 = src[i];
      int i2 = tgt[i];

      if (i1 < 0 || i1 >= data.nrow() || i2 < 0 || i2 >= data.nrow()) {
        gcc[i] = gcc_pval[i] = gcc1_vec[i] = gcc2_vec[i] =
          pcc[i] = pcc_pval[i] = scc[i] = scc_pval[i] =
          kcc[i] = kcc_pval[i] = NA_REAL;
        continue;
      }

      for (int j = 0; j < nSamples; j++) {
        xData[j] = data(i1, j);
        yData[j] = data(i2, j);
        xRank[j] = ranks(i1, j);
        yRank[j] = ranks(i2, j);
      }

      double gcc1_val = giniAsym(xData, yData);
      double gcc2_val = giniAsym(yData, xData);
      if (!R_finite(gcc1_val)) gcc1_val = 0.0;
      if (!R_finite(gcc2_val)) gcc2_val = 0.0;

      gcc1_vec[i] = gcc1_val;
      gcc2_vec[i] = gcc2_val;

      double chosen_gcc = asymmetricGCC ? gcc1_val : calcGiniSymmetric(xData, yData);
      if (!R_finite(chosen_gcc)) chosen_gcc = 0.0;
      gcc[i] = chosen_gcc;

      if (bootstrapIterations > 0) {
        gcc_pval[i] = permutationTests(chosen_gcc, xData, xRank, yData, yRank, wtVec, bootstrapIterations);
      } else {
        gcc_pval[i] = NA_REAL;
      }

      double r_pcc = pearsonCorrelation(xData, yData);
      pcc[i] = r_pcc;
      if (R_finite(r_pcc)) {
        double t_stat = r_pcc * sqrt((nSamples - 2) / (1 - r_pcc * r_pcc));
        pcc_pval[i] = 2 * R::pt(-fabs(t_stat), nSamples - 2, 1, 0);
      } else {
        pcc_pval[i] = NA_REAL;
      }

      double r_scc = spearmanCorrelation(xRank, yRank);
      scc[i] = r_scc;
      if (R_finite(r_scc)) {
        double t_stat = r_scc * sqrt((nSamples - 2) / (1 - r_scc * r_scc));
        scc_pval[i] = 2 * R::pt(-fabs(t_stat), nSamples - 2, 1, 0);
      } else {
        scc_pval[i] = NA_REAL;
      }

      double r_kcc = kendallCorrelation(xData, yData);
      kcc[i] = r_kcc;
      if (R_finite(r_kcc)) {
        double z_stat = r_kcc * sqrt((9 * nSamples * (nSamples - 1)) /
                                     (2 * (2 * nSamples + 5)));
        kcc_pval[i] = 2 * R::pnorm(-fabs(z_stat), 0, 1, 1, 0);
      } else {
        kcc_pval[i] = NA_REAL;
      }
    }
  }
};

// ===================== Exported function =====================
// [[Rcpp::export]]
DataFrame multi_corr_necorr(NumericMatrix expression,
                            IntegerMatrix ranks,
                            IntegerVector src,
                            IntegerVector tgt,
                            int bootstrapIterations = 0,
                            bool useBestGCC = false,
                            bool asymmetricGCC = false,
                            CharacterVector rownames_in = CharacterVector()) {
  int nEdges = src.size();
  int nSamples = expression.ncol();

  NumericVector wt(nSamples);
  for (int j = 0; j < nSamples; j++) {
    wt[j] = (2 * (j + 1) - nSamples - 1);
  }

  NumericVector gcc(nEdges), gcc_pval(nEdges),
  gcc1_vec(nEdges), gcc2_vec(nEdges),
  pcc(nEdges), pcc_pval(nEdges),
  scc(nEdges), scc_pval(nEdges),
  kcc(nEdges), kcc_pval(nEdges);

  MultiCorrWorker worker(expression, ranks, wt, src, tgt, nSamples,
                         bootstrapIterations, useBestGCC, asymmetricGCC,
                         gcc, gcc_pval, gcc1_vec, gcc2_vec,
                         pcc, pcc_pval, scc, scc_pval, kcc, kcc_pval);

  parallelFor(0, nEdges, worker);

  DataFrame df = DataFrame::create(
    _["GCC"] = gcc, _["GCC_pvalue"] = gcc_pval,
    _["GCC1"] = gcc1_vec, _["GCC2"] = gcc2_vec,
    _["PCC"] = pcc, _["PCC_pvalue"] = pcc_pval,
    _["SCC"] = scc, _["SCC_pvalue"] = scc_pval,
    _["KCC"] = kcc, _["KCC_pvalue"] = kcc_pval
  );

  if (rownames_in.size() == nEdges) {
    df.attr("row.names") = rownames_in;
  }

  return df;
}
