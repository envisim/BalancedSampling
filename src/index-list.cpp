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

IndexList* IndexList::copy() {
  IndexList* il = new IndexList(olen);

  for (int i = 0; i < olen; i++) {
    il->list[i] = list[i];
    il->reverse[i] = reverse[i];
  }

  il->len = len;

  return il;
}

IndexList* IndexList::copyLen() {
  IndexList* il = new IndexList(olen);

  for (int i = 0; i < len; i++) {
    il->list[i] = list[i];
    il->reverse[list[i]] = i;
  }

  il->len = len;

  return il;
}


int IndexList::length() {
  return len;
}

void IndexList::fill() {
  for (int i = 0; i < olen; i++) {
    list[i] = i;
    reverse[i] = i;
  }

  len = olen;
}

void IndexList::reset() {
  len = olen;
  return;
}

void IndexList::set(const int id) {
  if (id < 0 || id >= olen)
    throw std::range_error("(set) Inadmissible value of id");

  list[id] = id;
  reverse[id] = id;
}

int IndexList::get(const int k) {
  if (k < 0 || k >= len)
    throw std::range_error("(get) Inadmissible value of k");

  return list[k];
}

int IndexList::getK(const int id) {
  if (id < 0 || id >= olen)
    throw std::range_error("(getK) Inadmissable value of id");

  return reverse[id];
}

bool IndexList::exists(const int id) {
  return reverse[id] < len && reverse[id] >= 0;
}

int IndexList::draw() {
  // int k = (int)((double) len * stduniform());
  int k;
  do { k = (int)((double)(len) * stduniform()); } while (k == len);
  return list[k];
}

void IndexList::erase(const int id) {
  if (id < 0 || id >= olen)
    throw std::range_error("(erase, 1) Inadmissible value of id");

  int k = reverse[id];
  if (k < 0 || k >= len)
    throw std::range_error("(erase, 2) Inadmissible value of id, k");

  len -= 1;

  // Early return, no need to swap
  if (k == len)
    return;

  list[k] = list[len];
  list[len] = id;

  reverse[list[k]] = k;
  reverse[id] = len;
}

void IndexList::eraseK(const int k) {
  if (k < 0 || k >= len)
    throw std::range_error("(eraseK) Inadmissible value of k");

  len -= 1;

  // Early return, no need to swap
  if (k == len)
    return;

  int id = list[k];

  list[k] = list[len];
  list[len] = id;

  reverse[list[k]] = k;
  reverse[id] = len;
}
