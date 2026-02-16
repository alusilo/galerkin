//=============================
// Point source geometry
//=============================

x_max = 3000.; // max width
y_max = 3000.; // max heigth

mpl_layer = 300.;

// Characteristic length for
// * simulation area
// len <= vs/(3*freq), 3 is the number of triangles per wave length
len = 82.;

Point(1) = {0.0, 0.0, 0.0, len};
Point(2) = {0.0, y_max, 0.0, len};
Point(3) = {x_max, y_max, 0.0, len};
Point(4) = {x_max, 0.0, 0.0, len};
Point(5) = {x_max-mpl_layer, 0.0, 0.0, len};
Point(6) = {x_max-mpl_layer, y_max-mpl_layer, 0.0, len};
Point(7) = {mpl_layer, y_max-mpl_layer, 0.0, len};
Point(8) = {mpl_layer, 0.0, 0.0, len};

// Define the physical points
Physical Point(0) = {1,2,3,4,5,6,7,8};

// define the edges
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,1};
Line(9) = {8,5};

// Define the physical lines
Physical Line(0) = {1,2,3,4,5,6,7,8,9};

// define the surface
Line Loop(1) = {-5,-6,-7,-9};
Line Loop(2) = {1,2,3,4,5,6,7,8};

Plane Surface(1) = {1};
Plane Surface(2) = {2};

// Define the rest of the surfaces
Physical Surface(0) = {1,2};
