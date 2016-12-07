#include <Rcpp.h>
#include <algorithm>
#include <map>
#include <bitset>
#define NBITS 120

using namespace Rcpp;
using namespace std;

double calcGCC(vector<double> &xData, vector<int> &xIdx, vector<int> &yIdx, vector<double> &wt);
double calcPvalue(vector<double> &xData, vector<int> &xIdx, vector<int> &yIdx, vector<double> &wt,
                  double theAbsOfTheRealGCC, int perm);

vector<int> indexValues(vector<double> &values) {
  int nSamples = values.size();
  vector<double> sortedValues (nSamples);
  partial_sort_copy(values.begin(),values.end(),sortedValues.begin(),sortedValues.end());
  vector<int> valuesIndex (nSamples);
  // iterate over values, find in sortedValues, store position
  vector<double>::iterator it;
  for(int k=0;k<nSamples;k++) {
    it = find(sortedValues.begin(),sortedValues.end(),values[k]);
    int rank = it - sortedValues.begin();
    while (valuesIndex[rank] > 0) {
      rank++;
    }
    valuesIndex[rank] = k+1;
  }
  for(int rank=0;rank<nSamples;rank++) {
    valuesIndex[rank]--;
  }
  return valuesIndex;
}

vector<double> r_to_cpp(NumericVector values) {
  vector<double> vec(values.size());
  for(int i=0;i<values.size();i++) {
    vec[i] = values[i];
  }
  return vec;
}

// [[Rcpp::export]]
DataFrame gini(DataFrame edges, DataFrame expression, int bootstrapIterations, double statCutoff) {
  // edges has two vectors: source and target
  // each is a list of gene identifiers
  // The expression data frame has one numeric vector per gene identifier

  // sort the expression data for each gene
  // identify null values (treat 0 as null)
  int nGenes = expression.size();
  int nSamples = expression.nrows();
  vector<string> geneNames = expression.names();
  map<string,int> idToOffset;
  vector<vector<double> > data;
  vector<vector<int> > sorted;
  vector<bitset<NBITS> > mask;
  sorted.reserve(nGenes);
  data.reserve(nGenes);
  for (int i=0; i<nGenes; i++) {
    idToOffset[geneNames[i]] = i;
    vector<double> vec = r_to_cpp(expression[i]);
    sorted.push_back(indexValues(vec));
    data.push_back(vec);
    bitset<NBITS> nullValues;
    for (int j=0;j<vec.size();j++) {
      if (vec[j] == 0) {
        nullValues[j]=1;
      }
    }
    mask.push_back(nullValues);
  }

  // calculate weight vector in advance
  vector<double> wt (nSamples);
  for(int j=0;j<nSamples;j++) {
    wt[j] = (double) 2*(j+1) - nSamples - 1;
  }

  // iterate over the source - target pairs
  int nPairs = edges.nrows();
  CharacterVector source = edges["source"];
  CharacterVector target = edges["target"];
  NumericVector gini, pvalue;
  for(int i=0;i<nPairs; i++) {
    string s = as<string>(source[i]);
    string t = as<string>(target[i]);
    double gcc = 0.0;
    double pval = 1.0;
    map<string,int>::iterator source_it = idToOffset.find(s);
    if (source_it != idToOffset.end()) {
      map<string,int>::iterator target_it = idToOffset.find(t);
      if (target_it != idToOffset.end()) {
        // process this pair
        int source_offset = source_it->second;
        int target_offset = target_it->second;
        if (mask[source_offset].count() > 0 || mask[target_offset].count() > 0) {
          bitset<NBITS> nulls = mask[source_offset] | mask[target_offset];
          int nSamples_nonNull = nSamples - nulls.count();
          if (nSamples_nonNull >= 4) {
            vector<double> sourceData (nSamples_nonNull);
            vector<double> targetData (nSamples_nonNull);
            vector<double> wt2 (nSamples_nonNull);
            int i=0;
            for(int j=0;j<nSamples;j++) {
              if (!nulls.test(j)) {
                wt2[i] = (double) 2*(i+1) - nSamples_nonNull - 1;
                sourceData[i] = data[source_offset][j];
                targetData[i] = data[target_offset][j];
                i++;
              }
            }

            vector<int> sourceIdx = indexValues(sourceData);
            vector<int> targetIdx = indexValues(targetData);

            gcc = calcGCC(sourceData,sourceIdx,targetIdx,wt2);
            if (abs(gcc) >= statCutoff) {
              // calculate p-value
              pval = calcPvalue(sourceData,sourceIdx,targetIdx,wt2,gcc,bootstrapIterations);
            }
          }
        }
        else {
          gcc = calcGCC(data[source_offset],sorted[source_offset],sorted[target_offset],wt);
          if (abs(gcc) >= statCutoff) {
            // calculate p-value
            pval = calcPvalue(data[source_offset],sorted[source_offset],sorted[target_offset],wt,gcc,bootstrapIterations);
          }
        }
      }
    }
    // append gcc and pval to edges dataframe
    gini.push_back(gcc);
    pvalue.push_back(pval);
  }
  return DataFrame::create(_["source"]= source, _["target"]= target, _["gini"]=gini, _["pval"]=pvalue);
}

double calcGCC(vector<double> &xData, vector<int> &xIdx, vector<int> &yIdx, vector<double> &wt) {
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

double calcPvalue(vector<double> &xData, vector<int> &xIdx, vector<int> &yIdx, vector<double> &wt,
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
