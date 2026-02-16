#include <Triangle.h>

#define MAX(a,b) (a >= b ? a : b)
#define MIN(a,b) (a <= b ? a : b)

class Mesh : public Triangle {
private:
  void read_triangles(string infile);
  void update_nodes(int n, float x, float y, int mark);
  void read_node(string infile);
public:
  int size;
  Triangle *tri;
  string filename;
  Mesh() {};
  Mesh(string ele_file, string node_file);
  //Mesh& operator =(const Mesh& m);
  Triangle& operator [](int idx);
  const Triangle& operator [](int idx) const;
  friend ostream& operator << (ostream& stream, const Mesh &m) {
    int i;
    for(i=0;i<m.size;i++)
      stream << m.tri[i];
    return stream;
  }
  void load_vel(float *alpha, float *beta, float *rho, int nx, int nz, float dx, float dz, float freq);
  void evaluate_mesh(float freq);
  void statistics();
  ~Mesh();
};
