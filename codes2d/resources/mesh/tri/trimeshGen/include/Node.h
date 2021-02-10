# include <cstdlib>
# include <iostream>
# include <sstream>
# include <iomanip>
# include <fstream>
# include <vector>
# include <cmath>
# include <cstring>

using namespace std;

class Node {
public:
  int n;
  int mark;
  double x;
  double y;
  Node() {};
  Node(int n, double x = 0., double y = 0., int mark = 0) : n(n), x(x), y(y), mark(mark) {};
  double distance(Node obj);
};
