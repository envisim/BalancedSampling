#ifndef KDNODE_HEADER
#define KDNODE_HEADER


class KDNode {
  // REGULAR NODE
public:
  KDNode *parent = nullptr;
  KDNode *cleft = nullptr; // <=
  KDNode *cright = nullptr; // >
  int split = -1;
  double value = 0.0;

private:
  int terminal = 0; // 1 if terminal node, 0 if regular node
  int nunits = 0;
public:
  int *units = nullptr;

public:
  KDNode(KDNode*, const int);
  ~KDNode();
  void copy(KDNode*);
  void prune(const int);

  void setTerminal(const int);
  int isTerminal();
  KDNode* getSibling();
  void addUnits(const int*, const int);
  void removeUnit(const int);
  bool exists(const int);
  int getSize();
};

#endif
