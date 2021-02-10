#include <Node.h>

class Triangle : public Node {
public:
  int m;
  float vp;
  float vs;
  float rho;
  Node node1;
  Node node2;
  Node node3;
  Triangle() {};
  //Triangle(float vel) : vel(vel) {};
  Triangle(int m, int n1, int n2, int n3, float x1=0., float y1=0., float x2=0., float y2=0., float x3=0., float y3=0.) : m(m), node1(n1,x1,y1), node2(n2,x2,y2), node3(n3,x3,y3) {};
  friend ostream& operator << (ostream& stream, const Triangle& t) {
    stream << "Triangle: " << t.m << " --> vel(" << t.vp << ")" << endl;
    stream << "\tNode: " << t.node1.n << endl;
    stream << "\t\t" << t.node1.x << " " << t.node1.y << " -> " << t.node1.mark << endl;
    stream << "\tNode: " << t.node2.n << endl;
    stream << "\t\t" << t.node2.x << " " << t.node2.y << " -> " << t.node2.mark << endl;
    stream << "\tNode: " << t.node3.n << endl;
    stream << "\t\t" << t.node3.x << " " << t.node3.y << " -> " << t.node3.mark << endl;

    return stream;
  }
  float longest();
  float shortest();
  float perim();
  float area();
  float cond();
  float cx();
  float cy();
};
