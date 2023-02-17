#include "unorderedset-draw.h"

#define intuniform(N) ((int)((double)N * stduniform()))

int unorderedsetDraw(std::unordered_set<int> &set, int N) {
  int k = intuniform(N);

  // If k exists in the set, then we choose it
  if (set.count(k)) {
    return k;
  }

  // Add k to the set
  set.insert(k);
  std::unordered_set<int>::iterator it = set.find(k);
  it++;
  // Select the bucket after k (or the first bucket, if no bucket after k)
  int idx = it == set.end() ? *set.begin() : *it;
  // Remove k again
  set.erase(k);

  return idx;
}
