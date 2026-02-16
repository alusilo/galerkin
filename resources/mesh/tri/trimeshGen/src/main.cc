#include <Mesh.h>
#include <iostream>
#include <sstream>
#include <string>
#include <valarray>
#include <rsf.hh>
using namespace std;

template <typename T>
string to_string ( T Number )
{
    ostringstream ss;
    ss << Number;
    return ss.str();
}

int main(int argc, char* argv[])
{
    sf_init(argc,argv); // Initialize RSF
    iRSF par(0), vp("vp"), vs("vs"), rh("rho");
    int n1, n2, d1, d2, size, i, iter;
    float *alpha, *beta, *rho;
    float freq;
    string filename="geo_model", nodesfile, elefile;

    par.get("freq", freq);
    par.get("iter", iter);
    vp.get("n1", n1);
    vp.get("n2", n2);
    vp.get("d1", d1);
    vp.get("d2", d2);

    size = n1*n2;

    elefile = filename+"."+to_string(iter)+".ele";
    nodesfile = filename+"."+to_string(iter)+".node";

    alpha = new float[size];
    beta = new float[size];
    rho = new float[size];

    valarray<float> _vp(size);
    valarray<float> _vs(size);
    valarray<float> _rh(size);

    vp >> _vp;
    vs >> _vs;
    rh >> _rh;

    for(i=0;i<size;i++) {
        alpha[i] = _vp[i];
        beta[i]  = _vs[i];
        rho[i]   = _rh[i];
    }

    Mesh m(elefile, nodesfile);
    m.load_vel(alpha, beta, rho, n2, n1, d2, d1, freq);

    delete rho;
    delete beta;
    delete alpha;

    exit(0);
}
