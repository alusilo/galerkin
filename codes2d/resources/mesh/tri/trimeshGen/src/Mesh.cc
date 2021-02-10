#include <Mesh.h>

Mesh::Mesh(string ele_file, string node_file) {
  this->read_triangles(ele_file);
  this->read_node(node_file);
  this->filename = ele_file+".param";
}
Triangle& Mesh::operator [](int idx) { return tri[idx]; }
const Triangle& Mesh::operator [](int idx) const { return tri[idx]; }
Mesh::~Mesh() {
  delete tri;
}
void Mesh::read_triangles(string infile) {
  ifstream inunit;
  int t_number, n_number, attr_number, m, n1, n2, n3;

  inunit.open ( infile.c_str ( ) );

  if (inunit.good()) {
    int i;
    string line;
    getline(inunit,line);
    istringstream is(line);
    is >> t_number >> n_number >> attr_number;

    this->size = t_number;

    this->tri = new Triangle[t_number];

    for(i=0;i<t_number;i++) {
      getline(inunit,line);
      istringstream is(line);
      is >> m >> n1 >> n2 >> n3;
      this->tri[i] = Triangle(m, n1, n2, n3);
    }
  } else {
    cerr << infile << ": Bad input file!";
    exit(1);
  }

  inunit.close ( );
}
void Mesh::update_nodes(int n, float x, float y, int mark) {
  int i;
  for(i=0;i<this->size;i++) {
    if (this->tri[i].node1.n == n)
      this->tri[i].node1 = Node(n, x, y, mark);
    if (this->tri[i].node2.n == n)
      this->tri[i].node2 = Node(n, x, y, mark);
    if (this->tri[i].node3.n == n)
      this->tri[i].node3 = Node(n, x, y, mark);
  }
}
void Mesh::read_node(string infile) {
  ifstream inunit;
  float x, y;
  int n, n_number, dimension, attr_number, attr, markers, mark;
  int n1 = this->node1.n, n2 = this->node2.n, n3 = this->node3.n;

  inunit.open ( infile.c_str ( ) );

  if (inunit.good()) {
    int i;
    string line;
    getline(inunit,line);
    istringstream is(line);
    is >> n_number >> dimension >> attr_number >> markers;
    for(i=0;i<n_number;i++) {
      getline(inunit,line);
      istringstream is(line);
      is >> n >> x >> y >> mark;
      this->update_nodes(n,x,y,mark);
    }
  } else {
    cerr << infile << ": Bad input file!";
    exit(1);
  }

  inunit.close();
}
void Mesh::load_vel(float *alpha, float *beta, float *rho, int nx, int nz, float dx, float dz, float freq) {
  FILE *output = fopen(this->filename.c_str(), "w");
  float param, lop, ld, p1, p2;
  int i,j,ix,iz,idx;
  cout << this->size << endl;
  fprintf(output, "# density[kg/m^3] s-wave velocity [m/s] p-wave velocity [m/s]\n");
  for(i=0;i<this->size;i++){
    ix = (int)((this->tri[i].cx())/dx);
    iz = (int)((this->tri[i].cy())/dz);

    idx = ix*nz + iz;

    this->tri[i].vp  = alpha[idx];
    this->tri[i].vs  = beta[idx];
    this->tri[i].rho = rho[idx];

    param = this->tri[i].area()/this->tri[i].perim();
    ld = this->tri[i].vp/freq;
    lop = ld/30.;

    if(lop < param)
      cout << tri[i].m << " " << 0.9998*tri[i].area() << endl;
    else
      cout << tri[i].m << " " << -1 << endl;

    fprintf(output, "%f %f %f\n", this->tri[i].rho, this->tri[i].vs, this->tri[i].vp);
  }
  fclose(output);
}

void Mesh::evaluate_mesh(float freq) {

}

void Mesh::statistics() {
  int i;
  float min_length=this->tri[0].shortest(), max_length=this->tri[0].longest(), min_area = this->tri[0].area(), max_area = this->tri[0].area();
  for(i=1;i<this->size;i++) {
    min_length = MIN(min_length, this->tri[i].shortest());
    max_length = MAX(max_length, this->tri[i].longest());
    min_area = MIN(min_area, this->tri[i].area());
    max_area = MAX(max_area, this->tri[i].area());
  }
  cout << "====================== Mesh ========================\n";
  cout << "Longest trinagle side: " << max_length << endl;
  cout << "Shortest trinagle side: " << min_length << endl;
  cout << "Maximun trinagle area: " << max_area << endl;
  cout << "Minimun trinagle area: " << min_area << endl;
  cout << "====================================================\n";
}
