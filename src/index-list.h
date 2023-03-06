#ifndef INDEXLIST_HEADER
#define INDEXLIST_HEADER

#include "uniform.h"

class IndexList {
private:
  int *list = nullptr;
  int *reverse = nullptr;
  int len = 0;
  int olen = 0;
public:
  IndexList(const int);
  ~IndexList();
  int length();
  void fill();
  void set(const int);
  int get(const int);
  bool exists(const int);
  int draw();
  void erase(const int);
  void eraseK(const int);
};

#endif
