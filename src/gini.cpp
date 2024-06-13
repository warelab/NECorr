#include <Rcpp.h>
#include <RcppParallel.h>
#include <iostream>
#include <algorithm>
#include <random>
#include <map>

//[[Rcpp::depends(RcppParallel)]]
using namespace Rcpp;
using namespace RcppParallel;
using namespace std;

// worker to rank the expression data in parallel
struct expressionRank : public Worker {
  const RMatrix<double> data;
  RMatrix<int> rank;

  expressionRank(const NumericMatrix data, IntegerMatrix rank)
    : data(data), rank(rank) {
  }

  void operator()(size_t begin, size_t end) {
    for(size_t i = begin; i < end; i++) {
      RMatrix<double>::Row row = data.row(i);
      RMatrix<int>::Row rankRow = rank.row(i);

      vector<double> sortedRow(row.length());
      partial_sort_copy(row.begin(),row.end(),sortedRow.begin(),sortedRow.end());

      for(int k=0;k<row.length();k++) {
        vector<double>::iterator it = find(sortedRow.begin(),sortedRow.end(),row[k]);
        int r = it - sortedRow.begin();
        while (rankRow[r] > 0) {
          r++;
        }
        rankRow[r] = k+1;
      }
      for(int r=0; r<row.length(); r++) {
        rankRow[r]--;
      }
    }
  }
};

// worker to find pairs in parallel
struct expressionOffsets : public Worker {
  const CharacterVector sourceNames, targetNames;
  map<string,int> nameToOffset;
  IntegerVector source, target;

  expressionOffsets(const CharacterVector sourceNames, const CharacterVector targetNames,
                    map<string,int> nameToOffset,
                    IntegerVector source, IntegerVector target)
    : sourceNames(sourceNames), targetNames(targetNames), nameToOffset(nameToOffset), source(source), target(target)  {
  }

  void operator()(size_t begin, size_t end) {
    // iterate over the pairs and store the offset or -1

    for(size_t i = begin; i < end; i++) {
      string s = as<string>(sourceNames[i]);
      map<string,int>::iterator it = nameToOffset.find(s);
      source[i] = (it == nameToOffset.end()) ? -1 : it->second;
      string t = as<string>(targetNames[i]);
      it = nameToOffset.find(t);
      target[i] = (it == nameToOffset.end()) ? -1 : it->second;
    }
  }
};

double calcGCC(vector<double> &xData, vector<int> &xIdx, vector<int> &yIdx, vector<double> &wt, int n) {
  double numerator=0;
  double denominator=0;
  for(int i=0;i<n;i++) {
    numerator += wt[i] * xData[yIdx[i]];
    denominator += wt[i] * xData[xIdx[i]];
  }
  return denominator > 0 ? numerator/denominator : 0;
}
static unsigned int x=12345, y=65432, z=362436069, w=521288629;
unsigned int xorshift128(void) {
  unsigned int t = x;
  t ^= t << 11;
  t ^= t >> 8;
  x = y; y = z; z = w;
  w ^= w >> 19;
  w ^= t;
  return w;
}

std::random_device rd;
std::mt19937_64 rng(rd());

double permutationTest(double originalGCC, vector<double> &xData, vector<int> &xIdx, vector<int> &yIdx, vector<double> &wt, int perm) {
  int n = xData.size();
  vector<int> vec(yIdx);

  // precompute some parts of the GCC
  double denominator=0;
  for(int i=0; i<n; i++) {
    denominator += wt[i] * xData[xIdx[i]];
  }
  if (denominator == 0) {
    return 1.0;
  }
  int m = 0; // number of times subsample gcc is more extreme than originalGCC
  for(int p=0; p<perm; p++) {
    std::shuffle(vec.begin(), vec.end(), rng);
    double numerator=0;
    for(int i=0; i<n; i++) {
      numerator += wt[i] * xData[vec[i]];
    }
    double gcc = numerator/denominator;
    if (originalGCC > 0) {
      if (gcc > originalGCC) {
        m++;
      }
    }
    else if (gcc < originalGCC) {
      m++;
    }
    if ((p == 10 && m/p > 0.5)
          || (p == 100 && m/p > 0.2)
          || (p == 500 && m/p > 0.1)
       ) { // get out early
      perm = p;
    }
  }
  return (double)m/perm;
}

double confidenceIntervalAndPvalue(double originalGCC, vector<double> &xData, vector<int> &xIdx, vector<int> &yIdx, vector<double> &wt, int perm, double confidence, double &min, double &max) {
  vector<double> xData2 = xData;
  int n = xData.size();
  int m = 0; // number of times subsample gcc is more extreme than originalGCC
  vector<int> xRank (n);
  vector<int> yRank (n);
  vector<int> xIdx2 (n);
  vector<int> yIdx2 (n);

  vector<vector<int> > xRankPos (n);
  vector<vector<int> > yRankPos (n);

  vector<double> scoreBins2000 (2000);
  vector<double> scoreBins20 (20);

  for(int i=0;i<n;i++) {
    xRank[xIdx[i]] = i;
    yRank[yIdx[i]] = i;
  }
  for(int i=0; i<perm; i++) {
    for(int j=0; j < n; j++) {
      xRankPos[j].clear();
      yRankPos[j].clear();
    }
    for(int j=0; j < n; j++) {
      int k = xorshift128() % n; // choose a random offset
      xData2[j] = xData[k];
      xRankPos[xRank[k]].push_back(j);
      yRankPos[yRank[k]].push_back(j);
    }
    int xRankNew=0;
    int yRankNew=0;
    for(int j=0; j<n; j++) {
      if (xRankPos[j].size() > 0) {
        for(int k=0; k < xRankPos[j].size(); k++) {
          xIdx2[xRankNew++] = xRankPos[j][k];
        }
      }
      if (yRankPos[j].size() > 0) {
        for(int k=0; k < yRankPos[j].size(); k++) {
          yIdx2[yRankNew++] = yRankPos[j][k];
        }
      }
    }
    double gcc = calcGCC(xData2,xIdx2,yIdx2,wt,n);
    if (originalGCC > 0) {
      if (gcc > originalGCC) {
        m++;
      }
    }
    else if (gcc < originalGCC) {
      m++;
    }
    if (i > 10 && m/i > 0.3) { // get out early
      perm = i;
    }
    int bin = floor((gcc + 1.0)*1000);
    if (bin < 0) {
      bin = 0;
    }
    if (bin >= 2000) {
      bin = 1999;
    }
    scoreBins2000[bin]++;
    bin = floor((gcc + 1.0)*10);
    if (bin < 0) {
      bin = 0;
    }
    if (bin >= 20) {
      bin = 19;
    }
    scoreBins20[bin]++;
  }
  // instead of sorting 10000 numbers
  // find the first scoreBins20 that contains 2.5% of the data
  // then go to the scoreBins2000 below that and fine a more precise bin
  int lowerBound = floor(perm * (1.0 - confidence) / 2.0);
  int upperBound = perm - lowerBound;
  int bigBinStart=0;
  int bigBinTally=scoreBins20[0];
  while (bigBinTally < lowerBound) {
    bigBinStart++;
    bigBinTally += scoreBins20[bigBinStart];
  }
  int smallBinStart = 100 * bigBinStart;
  int smallBinTally = bigBinTally - scoreBins20[bigBinStart] + scoreBins2000[smallBinStart];
  while (smallBinTally < lowerBound) {
    smallBinStart++;
    smallBinTally += scoreBins2000[smallBinStart];
  }
  min = (double)smallBinStart/1000.0 - 1.0;
  int bigBinEnd = bigBinStart;
  while (bigBinTally < upperBound) {
    bigBinEnd++;
    bigBinTally += scoreBins20[bigBinEnd];
  }
  int smallBinEnd = 100 * bigBinEnd;
  smallBinTally = bigBinTally - scoreBins20[bigBinEnd] + scoreBins2000[smallBinEnd];
  while (smallBinTally < upperBound) {
    smallBinEnd++;
    smallBinTally += scoreBins2000[smallBinEnd];
  }
  max = (double)smallBinEnd/1000.0 - 1.0;
  return (double)m/perm;
}

// worker to calculate score edges in parallel
struct scoreEdges : public Worker {
  const IntegerVector source;
  const IntegerVector target;
  const RMatrix<double> data;
  const RMatrix<int> rank;
  const RVector<double> wt;
  const int bootstrapIterations;
  const double statCutoff;

  RVector<double> gini, pvalue, mingini, maxgini;

  scoreEdges(const IntegerVector source, const IntegerVector target,
             const NumericMatrix data, const IntegerMatrix rank, const NumericVector wt,
             const int bootstrapIterations, const double statCutoff,
             NumericVector gini, NumericVector pvalue, NumericVector mingini, NumericVector maxgini)
    : source(source), target(target), data(data), rank(rank), wt(wt),
      bootstrapIterations(bootstrapIterations), statCutoff(statCutoff),
      gini(gini), pvalue(pvalue), mingini(mingini), maxgini(maxgini) {}

  void operator()(size_t begin, size_t end) {
    for(size_t i=begin; i < end; i++) {
      gini[i] = 0.0;
      pvalue[i] = 1.0;
      mingini[i] = 0.0;
      maxgini[i] = 0.0;
      if (source[i] >= 0 && target[i] >= 0) {
        RMatrix<double>::Row sourceData = data.row(source[i]);
        RMatrix<double>::Row targetData = data.row(target[i]);
        RMatrix<int>::Row sourceRank = rank.row(source[i]);
        RMatrix<int>::Row targetRank = rank.row(target[i]);
        vector<double> sd, td, w;
        vector<int> sr, tr;
        int n = sourceData.length();
        for(size_t j=0;j<n;j++) {
          sd.push_back(sourceData[j]);
          td.push_back(targetData[j]);
          sr.push_back(sourceRank[j]);
          tr.push_back(targetRank[j]);
          w.push_back(wt[j]);
        }
        double gcc1 = calcGCC(sd, sr, tr, w, n);
        double gcc2 = calcGCC(td, tr, sr, w, n);
        if (abs(gcc1) > abs(gcc2)) {
          gini[i] = gcc1;
          if (abs(gcc1) >= statCutoff) {
            // pvalue[i] = confidenceIntervalAndPvalue(gini[i], sd, sr, tr, w, bootstrapIterations, 0.95, mingini[i], maxgini[i]);
            pvalue[i] = permutationTest(gini[i], sd, sr, tr, w, bootstrapIterations);
          }
        }
        else {
          gini[i] = gcc2;
          if (abs(gcc2) >= statCutoff) {
            // pvalue[i] = confidenceIntervalAndPvalue(gini[i], td, tr, sr, w, bootstrapIterations, 0.95, mingini[i], maxgini[i]);
            pvalue[i] = permutationTest(gini[i], td, tr, sr, w, bootstrapIterations);
          }
        }
      }
    }
  }
};

//' Do values in a numeric vector fall in specified range?
//'
//' This is a shortcut for `x >= left & x <= right`, implemented
//' efficiently in C++ for local values, and translated to the
//' appropriate SQL for remote tables.
//'
//' @param edges dataframe with the edges with 3 columns: source, interaction and
//' target
//' @param expression, expression table containing the genes
//' @param bootstrapIterations number of bootstrap
//' @param statCutoff stat cut-off p-value
//' @export
// [[Rcpp::export(rng = false)]]
DataFrame gini(DataFrame edges, NumericMatrix expression,
               int bootstrapIterations, double statCutoff) {
  // edges has (at least) two vectors: source and target
  // each is a list of gene identifiers
  // The expression data matrix has one row per gene, one column per sample

  int nSamples = expression.ncol();
  int nGenes = expression.nrow();
  int nPairs = edges.nrows();

  // calculate weight vector in advance
  NumericVector wt (nSamples);
  for(int j=0;j<nSamples;j++) {
    wt[j] = (double) 2*(j+1) - nSamples - 1;
  }

  // create a map for gene names
  List dimNames = expression.attr("dimnames");
  vector<string> geneNames = dimNames[0];
  map<string, int> nameToOffset;
  for(int i=0;i<nGenes;i++) {
    nameToOffset[geneNames[i]]=i;
  }

  // map gene source and target gene names to offsets in the expression matrix
  IntegerVector source(nPairs);
  IntegerVector target(nPairs);
  CharacterVector sourceNames = edges["source"];
  CharacterVector targetNames = edges["target"];
  expressionOffsets expressionOffsets(sourceNames, targetNames, nameToOffset, source, target);
  parallelFor(0,nPairs,expressionOffsets);

  // rank the samples for each gene
  IntegerMatrix rank(nGenes, nSamples);            // allocate space
  expressionRank expressionRank(expression, rank); // initialize worker
  parallelFor(0, nGenes, expressionRank);          // run with parallelFor

  // calculate gini correlation coefficient and p-value for each pair of genes
  NumericVector gini(nPairs);
  NumericVector pvalue(nPairs);
  NumericVector mingini(nPairs);
  NumericVector maxgini(nPairs);
  scoreEdges scoreEdges(source, target, expression, rank, wt, bootstrapIterations, statCutoff, gini, pvalue, mingini, maxgini);
  parallelFor(0, nPairs, scoreEdges);

  return DataFrame::create(_["source"]= sourceNames, _["target"]= targetNames, _["gini"]=gini, _["pvalue"]=pvalue, _["mingini"]=mingini, _["maxgini"]=maxgini);
}


