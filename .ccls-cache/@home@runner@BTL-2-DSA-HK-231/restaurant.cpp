#include "main.h"

static int MAXSIZE;

struct Node {
  int val;
  Node *left, *right;

  Node(int val, Node *left = NULL, Node *right = NULL) {
    this->val = val;
    this->left = left;
    this->right = right;
  }
};

class BSTree {
private:
  int size;

public:
  Node *root = nullptr;

  BSTree(const int val) {
    this->root = new Node(val);
    size = 1;
  }

  int getSize() { return this->size; }

  Node *insertNodeIterative(Node *root, const int val) {
    if (root == NULL) {
      root = new Node(val);
      return root;
    }
    Node *p = nullptr, *tmp = root;
    while (tmp) {
      p = tmp;
      if (val < tmp->val)
        tmp = tmp->left;
      else
        tmp = tmp->right;
    }
    if (val < p->val)
      p->left = new Node(val);
    else
      p->right = new Node(val);
    return root;
  }

  void insertNode(const int val) {
    insertNodeIterative(this->root, val);
    size++;
  }

  Node *delNodeRecur(Node *root, int val) {
    if (!root)
      return root;
    if (val < root->val) {
      root->left = delNodeRecur(root->left, val);
    } else if (val > root->val) {
      root->right = delNodeRecur(root->right, val);
    } else {
      if (!root->left && !root->right) {
        delete root;
        root = NULL;
        return root;
      } else if (!root->left) {
        Node *tmp = root->right;
        delete root;
        return tmp;
      } else if (!root->right) {
        Node *tmp = root->left;
        delete root;
        return tmp;
      } else {
        Node *tmp = root->right;
        while (tmp->left) {
          tmp = tmp->left;
        }
        root->val = tmp->val;
        root->right = delNodeRecur(root->right, tmp->val);
      }
    }
    return root;
  }

  void deleteNode(int val) {
    root = delNodeRecur(this->root, val);
    --size;
  }

  bool find(Node *root, int val) {
    if (!root)
      return false;
    if (root->val == val)
      return true;
    return find(root->left, val) || find(root->right, val);
  }
};

class MinHeap {
private:
  vector<int> heap;
  int maxsize;
  int count;

public:
  vector<vector<int>> value;
  vector<int> timest;

  MinHeap(int maxsize)
      : maxsize(maxsize), count(0) {}

  MinHeap(MinHeap *heap) {
    this->maxsize = heap->maxsize;
    this->count = heap->count;
    this->heap = heap->heap;
    this->value = heap->value;
    this->timest = heap->timest;
  }

  ~MinHeap() {}

  void resize(int maxsize) { this->maxsize = maxsize; }

  void insert(int val) {
    if (count == maxsize)
      return;
    count++;
    heap.push_back(val);
    reHeapUp(count - 1);
  }

  void insertFirst(int id, vector<int> res, int time) {
    if (count == maxsize)
      return;
    count++;
    heap.push_back(id);
    value.push_back(res);
    // for (int x : value[count - 1]) cout << x << " ";
    // cout << endl;
    timest.push_back(time);
    reHeapUp(count - 1);
  }

  bool isEmpty() { return !count; }

  int getSize() { return count; }

  int getPos(int id) {
    for (int i = 0; i < count; i++) {
      if (heap[i] == id)
        return i;
    }
    return -1;
  }

  int getID(int pos) { return heap[pos]; }

  bool find(int val) {
    return std::find(heap.begin(), heap.end(), val) != heap.end();
  }

  int top() { return heap[0]; }

  void pop() {
    if (count == 0)
      return;
    heap[0] = heap[count - 1];
    value[0] = value[count - 1];
    timest[0] = timest[count - 1];
    reHeapDown(0);
    count--;
  }

  void reHeapUp(int pos) {
    if (pos < 1)
      return;
    if (value[pos].size() < value[(pos - 1) / 2].size()) {
      swap(heap[pos], heap[(pos - 1) / 2]);
      value[pos].swap(value[(pos - 1) / 2]);
      swap(timest[pos], timest[(pos - 1) / 2]);
      reHeapUp((pos - 1) / 2);
    } else if (value[pos].size() == value[(pos - 1) / 2].size() &&
               timest[pos] < timest[(pos - 1) / 2]) {
      swap(heap[pos], heap[(pos - 1) / 2]);
      value[pos].swap(value[(pos - 1) / 2]);
      swap(timest[pos], timest[(pos - 1) / 2]);
      reHeapUp((pos - 1) / 2);
    }
  }
  void reHeapDown(int pos) {
    int leftChild = pos * 2 + 1, rightChild = pos * 2 + 2;
    int child = 0;
    if (leftChild < this->count) {
      if (rightChild < this->count &&
          value[rightChild].size() < value[leftChild].size()) {
        child = rightChild;
      } else if (rightChild < this->count &&
        value[rightChild].size() > value[leftChild].size()) {
        child = leftChild;
      }
      else if (rightChild < this->count &&
        value[rightChild].size() == value[leftChild].size()) {
        if (timest[rightChild] < timest[leftChild]) {
          child = rightChild;
        }
        else child = leftChild;
      }
      else child = leftChild;
    } else
      return;
    if (value[pos].size() > value[child].size()) {
      swap(heap[pos], heap[child]);
      value[pos].swap(value[child]);
      swap(timest[pos], timest[child]);
      reHeapDown(child);
    }
    else if (value[pos].size() == value[child].size() && timest[pos] > timest[child]) {
      swap(heap[pos], heap[child]);
      value[pos].swap(value[child]);
      swap(timest[pos], timest[child]);
      reHeapDown(child);
    }
  }

  void remove(int pos) {
    // TODO: remove the element with value equal to ID
    if (pos == count - 1) { 
      count--; 
      heap.pop_back();
      value.pop_back();
      timest.pop_back();
      return;
    }
    else {
      swap(heap[pos], heap[count - 1]);
      value[pos].swap(value[count - 1]);
      swap(timest[pos], timest[count - 1]);
      heap.pop_back();
      value.pop_back();
      timest.pop_back();
      count--;
      reHeapUp(pos);
      if (count != 0) reHeapDown(pos);
    }
  }

  void printHeap() {
    cout << "Heap: [";
    if (count == 0) {
      cout << "Empty" << endl;
      return;
    }
    for (int i = 0; i < count; ++i) {
      cout << heap[i] + 1 << " " << value[i].size() << " " << timest[i] << " - ";
    }
    cout << "]\n"; 
  }
};

class HuffmanNode {
private:
  char letter;
  int freq;
  bool isLeaf;
  int height;
  HuffmanNode *left, *right;
  friend class HuffmanTree;

public:
  HuffmanNode(char letter, int freq, bool isLeaf, int height,
              HuffmanNode *l = nullptr, HuffmanNode *r = nullptr)
      : letter(letter), freq(freq), isLeaf(isLeaf), height(height), left(l),
        right(r) {}
  virtual ~HuffmanNode(){};
  int getFreq() const { return this->freq; }
  char getletter() const { return this->letter; }
  bool checkLeaf() { return isLeaf; }
  int getHeight() { return this->height; }
  int getBalanceFac() {
    if (checkLeaf())
      return 0;
    return right->getHeight() - left->getHeight();
  }
  HuffmanNode *getLeftChild() { return this->left; }
  HuffmanNode *getRightChild() { return this->right; }
};

class HuffmanTree {
private:
public:
  HuffmanNode *root;
  int timestand;
  HuffmanTree(char letter, int freq, bool isLeaf, int timestand,
              HuffmanNode *l = nullptr, HuffmanNode *r = nullptr) {
    if (isLeaf) {
      this->root = new HuffmanNode(letter, freq, isLeaf, 1);
    } else {
      this->root = new HuffmanNode(
          letter, freq, isLeaf, 1 + max(l->getHeight(), r->getHeight()), l, r);
    }
    this->timestand = timestand;
  }

  void clearHuffNode(HuffmanNode *node) {
    if (node) {
      clearHuffNode(node->left);
      clearHuffNode(node->right);
    }
    delete node;
  }

  ~HuffmanTree() { clearHuffNode(root); }
  HuffmanNode *getRoot() { return this->root; }
  char getletter() { return root->getletter(); }
  int getFreq() { return root->getFreq(); }

  HuffmanNode *rotateLeft(HuffmanNode *curr) {
    if (curr == nullptr)
      return nullptr;
    HuffmanNode *subRight = curr->right;
    HuffmanNode *temp = subRight->left;
    subRight->left = curr;
    curr->right = temp;
    curr->height = max(curr->left->height, curr->right->height) + 1;
    subRight->height = max(subRight->left->height, subRight->right->height) + 1;
    return subRight;
  }

  HuffmanNode *rotateRight(HuffmanNode *curr) {
    if (curr == nullptr)
      return nullptr;
    HuffmanNode *subLeft = curr->left;
    HuffmanNode *temp = subLeft->right;
    subLeft->right = curr;
    curr->left = temp;
    curr->height = max(curr->left->height, curr->right->height) + 1;
    subLeft->height = max(subLeft->left->height, subLeft->right->height) + 1;
    return subLeft;
  }

  int getBalanceFac() {
    return root->right->getHeight() - root->left->getHeight();
  }

  bool unbalanceNodePreOrd(HuffmanNode *root) {
    if (!root)
      return false;
    if (root->getBalanceFac() < -1 || root->getBalanceFac() > 1)
      return true;
    return unbalanceNodePreOrd(root->left) || unbalanceNodePreOrd(root->right);
  }

  HuffmanNode *nodeUnbalancePreOrd(HuffmanNode *root) {
    if (!root)
      return root;
    if (root->getBalanceFac() < -1 || root->getBalanceFac() > 1)
      return root;
    if (!nodeUnbalancePreOrd(root->left)) {
      return nodeUnbalancePreOrd(root->right);
    }
    return nodeUnbalancePreOrd(root->left);
  }

  HuffmanNode *parent(HuffmanNode *root, HuffmanNode *child) {
    if (!root->left && !root->right)
      return nullptr;
    if (root->left == child || root->right == child) {
      return root;
    }
    if (!parent(root->left, child)) {
      return parent(root->right, child);
    }
    return parent(root->left, child);
  }

  void balanceTree(HuffmanNode *&root, int &step) {
    if (root == nullptr)
      return;
    while (unbalanceNodePreOrd(root) && step < 3) {
      HuffmanNode *x = nodeUnbalancePreOrd(root);
      HuffmanNode *y = parent(root, x);
      if (!y) {
        root = balanceNode(x, step);
      } else if (y->left == x) {
        y->left = balanceNode(x, step);
        y->height = 1 + max(y->left->height, y->right->height);
        while (parent(root, y)) {
          HuffmanNode *z = parent(root, y);
          z->height = 1 + max(z->left->height, z->right->height);
          y = z;
        }
      } else {
        y->right = balanceNode(x, step);
        y->height = 1 + max(y->left->height, y->right->height);
        while (parent(root, y)) {
          HuffmanNode *z = parent(root, y);
          z->height = 1 + max(z->left->height, z->right->height);
          y = z;
        }
      }
    }
  }

  HuffmanNode *balanceNode(HuffmanNode *root, int &step) {
    if (root->getBalanceFac() < -1 && root->left &&
        root->left->getBalanceFac() <= 0) {
      // LL
      root = rotateRight(root);
      step++;
    } else if (root->getBalanceFac() > 1 && root->right &&
               root->right->getBalanceFac() >= 0) {
      // RR
      root = rotateLeft(root);
      step++;
    } else if (root->getBalanceFac() < -1 && root->left &&
               root->left->getBalanceFac() > 0) {
      // LR
      root->left = rotateLeft(root->left);
      root = rotateRight(root);
      step++;
    } else if (root->getBalanceFac() > 1 && root->right &&
               root->right->getBalanceFac() < 0) {
      root->right = rotateRight(root->right);
      root = rotateLeft(root);
      step++;
    }
    return root;
  }
};

class Comparator {
public:
  bool operator()(HuffmanTree *h1, HuffmanTree *h2) const {
    if (h1->getFreq() != h2->getFreq()) {
      return h1->getFreq() < h2->getFreq();
    } else {
      if (h1->getletter() == ' ' && h2->getletter() != ' ') {
        return 0;
      } else if (h1->getletter() != ' ' && h2->getletter() == ' ') {
        return 1;
      } else if (h1->getletter() == ' ' && h2->getletter() == ' ') {
        return h1->timestand < h2->timestand;
      } else {
        return h1->timestand < h2->timestand;
      }
    }
  }
};

class Solution {
private:
  // HuffmanTree *hufftree;
  int sukunaTimeSt;
  MinHeap *sukunaRes;
  vector<BSTree *> gojoRes;
  vector<queue<int>> gojoResQ;
  vector<int> sukunaOrderId;
  stack<HuffmanTree *> ResSt; // luu thu tu LIFO cay huffman cua khach hang

public:
  Solution() {
    sukunaTimeSt = 0;
    sukunaRes = new MinHeap(MAXSIZE + 1);
    gojoRes.resize(MAXSIZE + 1);
    gojoResQ.resize(MAXSIZE + 1);
  }

  ~Solution() {
    delete sukunaRes;
    sukunaRes = nullptr;
    gojoRes.clear();
    gojoResQ.clear();
    sukunaOrderId.clear();
    // while (!ResSt.empty()) {
    //   delete ResSt.top();
    //   ResSt.pop();
    // }
  }

  // helper function for LAPSE method
  bool validName(string name, map<char, int> &mp) {
    for (char c : name)
      mp[c]++;
    return mp.size() >= 3;
  }

  static bool cmp1(pair<pair<char, int>, int> p1,
                   pair<pair<char, int>, int> p2) {
    if (p1.first.second == p2.first.second) {
      return p1.second < p2.second;
    }
    return p1.first.second < p2.first.second;
  }

  static bool cmp2(pair<pair<char, int>, int> p1,
                   pair<pair<char, int>, int> p2) {
    if (p1.first.second == p2.first.second) {
      if ((isupper(p1.first.first) && isupper(p2.first.first)) ||
          (islower(p1.first.first) && islower(p2.first.first)))
        return p1.first.first < p2.first.first;
      else
        return p1.first.first > p2.first.first;
    }
    return p1.first.second < p2.first.second;
  }

  char caesarCipher(char c, int num) {
    if (isupper(c)) {
      return (c - 'A' + num) % 26 + 'A';
    }
    return (c - 'a' + num) % 26 + 'a';
  }

  void caesarDecode(vector<pair<pair<char, int>, int>> &vp) {
    for (auto &x : vp) {
      char c = x.first.first;
      int shiftNum = x.first.second;
      char res = caesarCipher(c, shiftNum);
      x.first.first = res;
    }
  }

  void sortNameByFreq(vector<pair<pair<char, int>, int>> &vp, string name,
                      int n) {
    // Tinh tan suat xuat hien
    for (int i = 0; i < n; i++) {
      bool flag = false;
      int sz = vp.size();
      for (int j = 0; j < sz; j++) {
        if (vp[j].first.first == name[i]) {
          vp[j].first.second++;
          flag = true;
          break;
        }
      }
      if (!flag) {
        vp.push_back(make_pair(make_pair(name[i], 1), i));
      }
    }
    // sort danh sach theo thu tu tang dan tan suat xuat hien
    sort(vp.begin(), vp.end(), cmp1);
    // ma hoa danh sach
    caesarDecode(vp);
    // tien hanh cong don neu sau khi ma hoa co phan tu trung nhau
    int sz = vp.size();
    for (int i = sz - 1; i > 0; --i) {
      for (int j = i - 1; j >= 0; --j) {
        if (vp[i].first.first == vp[j].first.first) {
          vp[i].first.second += vp[j].first.second;
        }
      }
    }
    vector<bool> lowerchar(26, false);
    vector<bool> upperchar(26, false);
    vector<pair<pair<char, int>, int>> res;
    for (int i = sz - 1; i >= 0; --i) {
      if (!lowerchar[vp[i].first.first - 'a'] && islower(vp[i].first.first)) {
        lowerchar[vp[i].first.first - 'a'] = true;
        res.push_back(vp[i]);
      } else if (!upperchar[vp[i].first.first - 'A'] &&
                 isupper(vp[i].first.first)) {
        upperchar[vp[i].first.first - 'A'] = true;
        res.push_back(vp[i]);
      }
    }
    vp = res;
    for (auto x : vp) {
      cout << x.first.first << " " << x.first.second << endl;
    }
    // sort 1 lan nua, luu y theo thu tu alphabet va uu tien chu thuong truoc
    // chu hoa neu 2 ki tu co gia tri bang nhau
    sort(vp.begin(), vp.end(), cmp2);
    int sz1 = vp.size();
    for (int i = 0; i < sz1; ++i) {
      vp[i].second = i;
    }
  }

  HuffmanTree *buildHuffmanTree(vector<pair<pair<char, int>, int>> &vp, int n) {
    HuffmanTree *tmp1, *tmp2, *res = nullptr;
    set<HuffmanTree *, Comparator> pq;
    for (auto &x : vp) {
      HuffmanTree *tmp =
          new HuffmanTree(x.first.first, x.first.second, true, x.second);
      pq.insert(tmp);
    }
    if (pq.size() == 1) {
      return *pq.begin();
    }
    int timest = 0;
    while (pq.size() > 1) {
      tmp1 = *pq.begin();
      pq.erase(pq.begin());
      tmp2 = *pq.begin();
      pq.erase(pq.begin());
      res = new HuffmanTree(' ', tmp1->getFreq() + tmp2->getFreq(), false,
                            timest++, tmp1->getRoot(), tmp2->getRoot());
      int step = 0;
      res->balanceTree(res->root, step);
      // printHuffmanTreeTest(res->root, 0);
      // cout << endl;
      pq.insert(res);
    }
    // delete tmp1;
    // delete tmp2;
    return res;
  }

  void convertCharToBin(HuffmanNode *root, char c, string tmp, string &res) {
    if (!root)
      return;
    if (root->getletter() == c) {
      res = tmp;
      return;
    }
    convertCharToBin(root->getLeftChild(), c, tmp + "0", res);
    convertCharToBin(root->getRightChild(), c, tmp + "1", res);
  }

  void insertCusToGojoRes(int res, int id) {
    if (!gojoRes[id]) {
      gojoRes[id] = new BSTree(res);
      gojoResQ[id].push(res);
      return;
    }
    if (gojoRes[id] && !gojoRes[id]->root) {
      gojoRes[id] = new BSTree(res);
      gojoResQ[id].push(res);
      return;
    }
    gojoRes[id]->insertNode(res);
    gojoResQ[id].push(res);
  }

  void insertCusToSukunaRes(int res, int id) {
    id--;
    if (!sukunaRes->find(id)) {
      vector<int> v = {res};
      sukunaRes->insertFirst(id, v, sukunaTimeSt);
      sukunaOrderId.push_back(id);
      // sukunaRes->printHeap();
      // cout << endl;
      return;
    }
    int pos = sukunaRes->getPos(id);
    sukunaRes->value[pos].push_back(res);
    sukunaRes->timest[pos] = sukunaTimeSt;
    sukunaRes->reHeapDown(pos);
    if (sukunaOrderId.empty())
      sukunaOrderId.push_back(id);
    else {
      auto it = find(sukunaOrderId.begin(), sukunaOrderId.end(), id);
      if (it != sukunaOrderId.end()) {
        sukunaOrderId.erase(it);
        sukunaOrderId.push_back(id);
      } else {
        sukunaOrderId.push_back(id);
      }
    }
    // sukunaRes->printHeap();
    // cout << endl;
  }

  void LAPSE(string name) {
    int n = name.size();
    map<char, int> mp;
    if (!validName(name, mp))
      return;
    vector<pair<pair<char, int>, int>> vp;
    // Xay dung danh sach
    sortNameByFreq(vp, name, n);
    if (vp.size() == 0)
      return;
    // Xay dung cay Huffman tu danh sach
    if (vp.size() == 1) {
      insertCusToSukunaRes(0, 1);
      return;
    }
    for (auto x : vp) {
      cout << x.first.first << " " << x.first.second << endl;
    }
    HuffmanTree *hufftree = buildHuffmanTree(vp, n);
    ResSt.push(hufftree);
    // Dua vao cay Huffman da xay dung, chuyen doi cac ki tu sang ma nhi phan
    vector<pair<char, string>> huffmanCode;
    string name2 = "";
    for (int i = 0; i < n; i++) {
      int shiftnum = mp[name[i]];
      name2 += caesarCipher(name[i], shiftnum);
    }
    cout << name2 << endl;
    string binCode = "";
    for (int i = 0; i < n; i++) {
      string code = "", tmp = "";
      convertCharToBin(hufftree->getRoot(), name2[i], tmp, code);
      binCode += code;
    }
    // Neu binCode it hon hoac bang 10 ki tu thi thoi lay luon, con khong thi
    // lay 10 ki tu dau tien
    int m = binCode.size();
    cout << binCode << endl;
    if (m > 10) {
      binCode = binCode.substr(binCode.size() - 10, 10);
      m = 10;
    }
    int result = 0;
    // Chuyen day binCode sang so thap phan
    for (int i = m - 1; i >= 0; --i) {
      if (binCode[i] == '1')
        result += (1 << (i));
    }
    cout << result << endl;
    int id = result % MAXSIZE + 1;
    if (result & 1) {
      insertCusToGojoRes(result, id);
    } else {
      sukunaTimeSt++;
      insertCusToSukunaRes(result, id);
    }
  }

  // Helper function for KOKUSEN method
  void PostOrderArr(Node *root, vector<int> &v) {
    if (!root)
      return;
    PostOrderArr(root->left, v);
    PostOrderArr(root->right, v);
    v.push_back(root->val);
  }

  vector<vector<long long>> binCoeff;
  long long mod = 1e9 + 7;

  long long dfs(Node* root) {
    if (!root) return 1;

    int m = countNode(root), mleft = countNode(root->left);

    long long leftWays = dfs(root->left) % mod;
    long long rightWays = dfs(root->right) % mod;

    return (((leftWays * rightWays) % mod) *
            binCoeff[m - 1][mleft]) %
           mod;
  }

  int countNode(Node* root) {
    if (!root) return 0;
    return 1 + countNode(root->left) + countNode(root->right);
  }

  int countOfWays(Node* root) {
    int m = countNode(root);
    if (m == 1)
      return 1;
    // binCoeff of Pascal's triangle
    binCoeff.resize(m + 1);
    for (int i = 0; i < m + 1; ++i) {
      binCoeff[i] = vector<long long>(i + 1, 1);
      for (int j = 1; j < i; ++j) {
        binCoeff[i][j] = (binCoeff[i - 1][j - 1] + binCoeff[i - 1][j]) % mod;
      }
    }

    return ((dfs(root)) % mod) % MAXSIZE;
  }

  void CountAndDelBST(BSTree *tree, int i) {
    int step = countOfWays(tree->root);
    if (step > tree->getSize())
      step = tree->getSize();
    while (step--) {
      int x = gojoResQ[i].front();
      gojoResQ[i].pop();
      tree->deleteNode(x);
    }
  }

  void KOKUSEN() {
    int n = gojoRes.size();
    for (int i = 0; i < n; ++i) {
      if (gojoRes[i]) {
        CountAndDelBST(gojoRes[i], i);
      }
    }
  }

  // helper function for KEITEIKEN method
  void findIDToDel(int num, queue<pair<int, int>> &q) {
    int sz = sukunaRes->getSize();
    if (sz == 0)
      return;
    sz = min(num, sz);
    MinHeap *tmp = new MinHeap(sukunaRes);
    vector<int> ans;
    int i = 0;
    while (i < num && i < sz) {
      int top = tmp->top();
      tmp->pop();
      ans.push_back(top);
      i++;
    }
    for (int x : ans) {
      removeCus(x, num, q);
    }
  }

  void removeCus(int id, int num, queue<pair<int, int>> &q) {
    int pos = sukunaRes->getPos(id);
    if (num > sukunaRes->value[pos].size()) {
      num = sukunaRes->value[pos].size();
    }
    while (num-- && !sukunaRes->value[pos].empty()) {
      q.push(make_pair(*(sukunaRes->value[pos].begin()), id));
      sukunaRes->value[pos].erase(sukunaRes->value[pos].begin());
    }
    sukunaTimeSt++;
    if (sukunaRes->value[pos].empty()) {
      sukunaRes->remove(pos);
      return;
    }
    sukunaRes->timest[pos] = sukunaTimeSt;
    sukunaRes->reHeapUp(pos);
  }

  void KEITEIKEN(int num) {
    if (num == 0)
      return;
    queue<pair<int, int>> q;
    findIDToDel(num, q);
    while (!q.empty()) {
      auto top = q.front();
      cout << top.first << "-" << top.second + 1 << endl;
      q.pop();
    }
  }

  // helper function for HAND method
  void printHuffmanTree(HuffmanNode *root) {
    if (!root)
      return;
    printHuffmanTree(root->getLeftChild());
    if (root->getletter() == ' ') {
      cout << root->getFreq() << endl;
    } else {
      cout << root->getletter() << endl;
    }
    printHuffmanTree(root->getRightChild());
  }

  void HAND() {
    if (ResSt.empty())
      return;
    HuffmanNode *LastRoot = ResSt.top()->getRoot();
    printHuffmanTree(LastRoot);
    // delete LastRoot;
  }

  // helper function for LIMITLESS method
  void printInorder(Node *root) {
    if (!root)
      return;
    printInorder(root->left);
    cout << root->val << endl;
    printInorder(root->right);
  }

  void LIMITLESS(int num) {
    if (!gojoRes[num] || num > MAXSIZE)
      return;
    printInorder(gojoRes[num]->root);
  }

  // helper function for CLEAVE method
  void printPreorder(int pos, int num, int size) {
    if (pos >= size)
      return;
    int id = sukunaRes->getID(pos);
    int n = sukunaRes->value[pos].size();
    for (int i = n - 1, j = 0; i >= 0 && j < num; --i, ++j) {
      cout << id + 1 << "-" << sukunaRes->value[pos][i] << endl;
    }
    printPreorder(pos * 2 + 1, num, size);
    printPreorder(pos * 2 + 2, num, size);
  }

  void CLEAVE(int num) {
    int size = sukunaRes->getSize();
    printPreorder(0, num, size);
  }
};

// void simulate(string filename) {
//   ifstream ss(filename);
//   string str, name;
//   int maxsize, num;
//   ss >> str;
//   if (str == "MAXSIZE") {
//     ss >> maxsize;
//     MAXSIZE = maxsize;
//   }
//   Solution s;
//   int i = 2;
//   while (ss >> str) {
//     if (str == "LAPSE") // LAPSE <NAME>
//     {
//       ss >> name;
//       // cout << "Line " << i << endl;
//       i++;
//       // cout << "LAPSE " << endl;
//       s.LAPSE(name);
//     } else if (str == "KOKUSEN") // KOKUSEN
//     {
//       // cout << "Line " << i << endl;
//       i++;
//       // cout << "KOKUSEN " << endl;
//       s.KOKUSEN();
//     } else if (str == "KEITEIKEN") // KEITEIKEN <NUM>
//     {
//       // cout << "Line " << i << endl;
//       i++;
//       ss >> num;
//       // cout << "KEITEIKEN " << num << endl;
//       s.KEITEIKEN(num);
//     } else if (str == "HAND") // HAND
//     {
//       // cout << "Line " << i << endl;
//       i++;
//       // cout << "HAND " << endl;
//       s.HAND();
//     } else if (str == "LIMITLESS") // LIMITLESS <NUM>
//     {
//       // cout << "Line " << i << endl;
//       i++;
//       ss >> num;
//       // cout << "LIMITLESS " << num << endl;
//       s.LIMITLESS(num);
//     } else if (str == "CLEAVE") // DOMAIN_EXPANSION
//     {
//       // cout << "Line " << i << endl;
//       i++;
//       ss >> num;
//       // cout << "CLEAVE " << num << endl;
//       s.CLEAVE(num);
//     }
//   }
//   return;
// }


