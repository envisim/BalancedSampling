#include <algorithm>
#include <unordered_map>
#include <Rcpp.h>
#include "cube-internal.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

#define pbig(p, eps) ((p) >= 1.0 - (eps))
#define pclose(p, eps) ((p) <= (eps) || (p) >= 1.0 - (eps))
#define mati(r, c, p) ((r) * (p) + (c))

// [[Rcpp::export(.cube_stratified_cpp)]]
Rcpp::IntegerVector cube_stratified_cpp(
  Rcpp::NumericVector &prob,
  Rcpp::NumericMatrix &x,
  Rcpp::IntegerVector &strata,
  double eps
) {
  int N = x.nrow();
  int p = x.ncol();

  if (N != prob.length())
    std::range_error("prob does not match x");
  if (N != strata.length())
    std::range_error("strata does not match x");

  std::unordered_map<int, int> stratumMap;

  int sampleSize = 0;

  double *probabilities = new double[N];
  IndexList *idx = new IndexList(N);
  int *sample = new int[N];

  for (int i = 0; i < N; i++) {
    idx->set(i);
    probabilities[i] = prob[i];
  }

  for (int i = 0; i < N; i++) {
    if (pclose(prob[i], eps)) {
      if (pbig(prob[i], eps)) {
        sample[sampleSize] = i;
        sampleSize += 1;
      }

      idx->erase(i);

      continue;
    }

    if (stratumMap.count(strata[i]) == 0)
      stratumMap[strata[i]] = 1;
    else
      stratumMap[strata[i]] += 1;
  }

  idx->shuffle();

  int *tindex = new int[N];
  double *tprobabilities = new double[N];
  int *tsample = new int[N];
  int *tivec = new int[N];
  double *tuvec = new double[N];
  int *tstratumArr = new int[stratumMap.size()];
  int tstratumArrLen = 0;
  double *tamat = new double[N * (p + stratumMap.size())];
  double *tbmat = new double[(p + 1 + stratumMap.size()) * (p + stratumMap.size())];

  // 1. Flight per stratum
  int maxSize1 = p + 2;

  for (std::unordered_map<int, int>::iterator it = stratumMap.begin(); it != stratumMap.end(); ++it) {
    // We can skip this completely if there are too few obs in the stratum
    if (it->second < maxSize1)
      continue;

    int tsampleSize = 0;
    int stratLen = 0;
    IndexList* tidx = new IndexList(it->second);

    for (int i = 0; i < idx->length(); i++) {
      int id = idx->get(i);
      if (it->first != strata[id])
        continue;

      tidx->set(stratLen);
      tindex[stratLen] = id;
      tprobabilities[stratLen] = probabilities[id];
      tamat[mati(0, stratLen, it->second)] = 1.0;
      for (int k = 0; k < p; k++)
        tamat[mati(k + 1, stratLen, it->second)] = x(id, k) / probabilities[id];

      stratLen += 1;
    }

    cubeFlightPhase(
      tamat,
      tprobabilities,
      tidx,
      nullptr,
      tsample,
      &tsampleSize,
      eps,
      p + 1,
      it->second,
      tivec,
      nullptr,
      tuvec,
      tbmat
    );

    for (int i = 0; i < stratLen; i++) {
      int id = tindex[i];
      if (!tidx->exists(i)) {
        idx->erase(id);
        it->second -= 1;
      } else {
        probabilities[id] = tprobabilities[i];
      }
    }

    for (int i = 0; i < tsampleSize; i++) {
      sample[sampleSize] = tindex[tsample[i] - 1] + 1;
      sampleSize += 1;
    }

    if (it->second == 0) {
      tstratumArr[tstratumArrLen] = it->first;
      tstratumArrLen += 1;
    }

    delete tidx;
  }

  for (int k = 0; k < tstratumArrLen; k++) {
    stratumMap.erase(tstratumArr[k]);
  }

  tstratumArrLen = 0;
  for (std::unordered_map<int, int>::iterator it = stratumMap.begin(); it != stratumMap.end(); ++it) {
    tstratumArr[tstratumArrLen] = it->first;
    tstratumArrLen += 1;
  }

  if (tstratumArrLen != stratumMap.size())
    std::range_error("stratum map error");

  // 2. Flight on remaining full population
  int maxSize2 = p + 1 + stratumMap.size();

  if (idx->length() >= maxSize2) {
    int idxlen = idx->length();

    for (int i = 0; i < idxlen; i++) {
      int id = idx->get(i);
      tindex[i] = id;

      for (int k = 0; k < stratumMap.size(); k++)
        tamat[mati(k, id, N)] = strata[id] == tstratumArr[k] ? 1.0 : 0.0;

      for (int k = 0; k < p; k++)
        tamat[mati(k + stratumMap.size(), id, N)] = x(id, k) / prob[id];
    }

    cubeFlightPhase(
      tamat,
      probabilities,
      idx,
      nullptr,
      sample,
      &sampleSize,
      eps,
      p + stratumMap.size(),
      N,
      tivec,
      nullptr,
      tuvec,
      tbmat
    );

    for (int i = 0; i < idxlen; i++) {
      int id = tindex[i];
      if (idx->exists(id))
        continue;
      stratumMap[strata[id]] -= 1;
    }

    tstratumArrLen = 0;
    for (std::unordered_map<int, int>::iterator it = stratumMap.begin(); it != stratumMap.end(); ++it) {
      if (it->second > 0)
        continue;
      else if (it->second < 0)
        std::range_error("too many removed");

      tstratumArr[tstratumArrLen] = it->first;
      tstratumArrLen += 1;
    }

    for (int k = 0; k < tstratumArrLen; k++) {
      stratumMap.erase(tstratumArr[k]);
    }
  }

  // 3. Landing per stratum
  // int maxSize3 = p + 2;

  for (std::unordered_map<int, int>::iterator it = stratumMap.begin(); it != stratumMap.end(); ++it) {
    int tsampleSize = 0;
    int stratLen = 0;
    IndexList* tidx = new IndexList(it->second);

    for (int i = 0; i < idx->length(); i++) {
      int id = idx->get(i);
      if (it->first != strata[id])
        continue;

      tidx->set(stratLen);
      tindex[stratLen] = id;
      tprobabilities[stratLen] = probabilities[id];
      tamat[mati(0, stratLen, it->second)] = 1.0;
      for (int k = 0; k < p; k++)
        tamat[mati(k + 1, stratLen, it->second)] = x(id, k) / probabilities[id];

      stratLen += 1;
    }

    cubeLandingPhase(
      tamat,
      tprobabilities,
      tidx,
      // nullptr,
      tsample,
      &tsampleSize,
      eps,
      p + 1,
      it->second,
      tivec,
      tuvec,
      tbmat
    );

    for (int i = 0; i < tsampleSize; i++) {
      sample[sampleSize] = tindex[tsample[i] - 1] + 1;
      sampleSize += 1;
    }

    delete tidx;
  }

  // Clean up
  std::sort(sample, sample + sampleSize);
  Rcpp::IntegerVector svec(sample, sample + sampleSize);

  delete[] probabilities;
  delete idx;
  delete[] sample;

  delete[] tindex;
  delete[] tprobabilities;
  delete[] tsample;
  delete[] tivec;
  delete[] tuvec;
  delete[] tstratumArr;
  delete[] tamat;
  delete[] tbmat;

  return svec;
}