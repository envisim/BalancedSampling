#include "index-list.h"

IndexList::IndexList(const int N) {
  list = new int[N];
  reverse = new int[N];
  len = N;
  olen = N;
}

IndexList::~IndexList() {
  delete[] list;
  delete[] reverse;
}

int IndexList::length() {
  return len;
}

void IndexList::fill() {
  for (int i = 0; i < len; i++) {
    list[i] = i;
    reverse[i] = i;
  }
}

void IndexList::set(const int id) {
  list[id] = id;
  reverse[id] = id;
}

int IndexList::get(const int k) {
  if (k < 0 || k >= len)
    throw std::range_error("Inadmissible value of k");

  return list[k];
}

bool IndexList::exists(const int id) {
  return reverse[id] < len && reverse[id] >= 0;
}

int IndexList::draw() {
  int k = (int)((double) len * stduniform());
  return list[k];
}

void IndexList::erase(const int id) {
  if (id < 0 || id >= olen)
    throw std::range_error("Inadmissible value of id");

  int k = reverse[id];
  if (k < 0 || k >= len)
    throw std::range_error("Inadmissible value of id");

  len -= 1;

  list[k] = list[len];
  list[len] = id;

  reverse[list[k]] = k;
  reverse[id] = len;
}

void IndexList::eraseK(const int k) {
  if (k < 0 || k >= len)
    throw std::range_error("Inadmissible value of k");

  len -= 1;

  int id = list[k];

  list[k] = list[len];
  list[len] = id;

  reverse[list[k]] = k;
  reverse[id] = len;
}
