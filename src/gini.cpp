// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <algorithm>
#include <map>
#include <bitset>
#define NBITS 120

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

inline double calcGCC(vector<double> &xData, vector<int> &xIdx, vector<int> &yIdx, vector<double> &wt) {
  double numerator=0;
  double denominator=0;
  for(int i=0;i<xData.size();i++) {
    numerator += wt[i] * xData[yIdx[i]];
    denominator += wt[i] * xData[xIdx[i]];
  }
  return numerator/denominator;
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

void shuffle(vector<double> &data, vector<int> &idx, vector<int> &rank) {
  int n = data.size();
  for(int i=n-1;i>0;i--) {
    // int j = rand() % i;
    int j = xorshift128() % i;
    int tmpD = data[i]; data[i] = data[j]; data[j] = tmpD;
    int tmpI = rank[i]; rank[i] = rank[j]; rank[j] = tmpI;
  }
  for(int i=0;i<n;i++) {
    idx[rank[i]] = i;
  }
}

inline double calcPvalue(vector<double> &xData, vector<int> &xIdx, vector<int> &yIdx, vector<double> &wt,
                  double theRealGCC, int perm) {
  // make copies of xData and xIdx
  vector<double> rData = xData;
  vector<int> rIdx = xIdx;
  int n = rData.size();
  // make a rank vector
  vector<int> rRank (n);
  for(int i=0;i<n;i++) {
    rRank[rIdx[i]] = i;
  }
  // x = rand();
  // y = rand();
  // z = rand();
  // w = rand();
  int m=0;
  for(int i=1;i<=perm;i++) {
    shuffle(rData,rIdx,rRank);
    double gcc = calcGCC(rData,rIdx,yIdx,wt);
    if (theRealGCC > 0) {
      if (gcc > theRealGCC) {
        m++;
        if (i < 10 && m > 3) {
          return (double) m/i;
        }
        else if (i < 100 && m > 20) {
          return (double) m/i;
        }
      }
    }
    else {
      if (gcc < theRealGCC) {
        m++;
        if (i < 10 && m > 3) {
          return (double) m/i;
        }
        else if (i < 100 && m > 20) {
          return (double) m/i;
        }
      }
    }
  }
  if (m==0) {
    m = 1;
  }
  return (double) m/perm;
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

  RVector<double> gini, pvalue;

  scoreEdges(const IntegerVector source, const IntegerVector target,
             const NumericMatrix data, const IntegerMatrix rank, const NumericVector wt,
             const int bootstrapIterations, const double statCutoff,
             NumericVector gini, NumericVector pvalue)
    : source(source), target(target), data(data), rank(rank), wt(wt),
      bootstrapIterations(bootstrapIterations), statCutoff(statCutoff),
      gini(gini), pvalue(pvalue) {}

  void operator()(size_t begin, size_t end) {
    for(size_t i=begin; i < end; i++) {
      double gcc = 0.0;
      double pval = 1.0;
      if (source[i] >= 0 && target[i] >= 0) {
        RMatrix<double>::Row sourceData = data.row(source[i]);
        RMatrix<int>::Row sourceRank = rank.row(source[i]);
        RMatrix<int>::Row targetRank = rank.row(target[i]);
        vector<double> sd, w;
        vector<int> sr, tr;
        for(size_t j=0;j<sourceData.length();j++) {
          sd.push_back(sourceData[j]);
          sr.push_back(sourceRank[j]);
          tr.push_back(targetRank[j]);
          w.push_back(wt[j]);
        }
        gcc = calcGCC(sd, sr, tr, w);
        if (abs(gcc) >= statCutoff) {
          pval = calcPvalue(sd, sr, tr, w, gcc, bootstrapIterations);
        }
      }
      gini[i] = gcc;
      pvalue[i] = pval;
    }
  }
};

// [[Rcpp::export]]
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

  // map gene source and target gene names to offsets in the expression matrix
  IntegerVector source(nPairs);
  IntegerVector target(nPairs);

  List dimNames = expression.attr("dimnames");
  vector<string> geneNames = dimNames[0];
  map<string, int> nameToOffset;
  for(int i=0;i<nGenes;i++) {
    nameToOffset[geneNames[i]]=i;
  }

  // rank the samples for each gene
  IntegerMatrix rank(nGenes, nSamples);            // allocate space
  expressionRank expressionRank(expression, rank); // initialize worker
  parallelFor(0, nGenes, expressionRank);          // run with parallelFor

  // iterate over the pairs and store the offset or -1
  CharacterVector sourceNames = edges["source"];
  CharacterVector targetNames = edges["target"];
  map<string,int>::iterator it;
  for(int i=0;i<nPairs;i++) {
    string s = as<string>(sourceNames[i]);
    it = nameToOffset.find(s);
    source[i] = (it == nameToOffset.end()) ? -1 : it->second;
    string t = as<string>(targetNames[i]);
    it = nameToOffset.find(t);
    target[i] = (it == nameToOffset.end()) ? -1 : it->second;
  }

  // calculate gini correlation coefficient and p-value for each pair of genes
  NumericVector gini(nPairs);
  NumericVector pvalue(nPairs);
  scoreEdges scoreEdges(source, target, expression, rank, wt, bootstrapIterations, statCutoff, gini, pvalue);
  parallelFor(0, nPairs, scoreEdges);

  return DataFrame::create(_["source"]= sourceNames, _["target"]= targetNames, _["gini"]=gini, _["pval"]=pvalue);
}


