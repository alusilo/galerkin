//=============================
// Point source geometry
//=============================

x_max = 7600.; // max width
y_max = 5600.; // max heigth

pml_layer = 300.0;

// Characteristic length for
// * simulation area
// len <= vs/(3*freq), 3 is the number of triangles per wave length
vp = 4000.0;
vs = 2310.0;//vp/1.8;
freq = 7.5;
len = vs/(3*freq);

Point(1) = {0.0, 0.0, 0.0, len};
Point(2) = {pml_layer, 0.0, 0.0, len};
Point(3) = {x_max-pml_layer, 0.0, 0.0, len};
Point(4) = {x_max, 0.0, 0.0, len};
Point(5) = {x_max, y_max, 0.0, len};
Point(6) = {0.0, y_max, 0.0, len};
Point(7) = {pml_layer, pml_layer, 0.0, len};
Point(8) = {x_max-pml_layer, pml_layer, 0.0, len};
Point(9) = {x_max-pml_layer, y_max-pml_layer, 0.0, len};
Point(10) = {pml_layer, y_max-pml_layer, 0.0, len};

// Define the physical points
Physical Point(0) = {1,2,3,4,5,6,7,8,9,10};

// define the edges
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};
Line(7) = {2,7};
Line(8) = {3,8};
Line(9) = {7,8};
Line(10) = {8,9};
Line(11) = {9,10};
Line(12) = {10,7};

// Define the physical lines
Physical Line(0) = {1,2,3,4,5,6,7,8,9,10,11,12};

// define the surface
Line Loop(1) = {7,-12,-11,-10,-8,3,4,5,6,1};
Line Loop(2) = {9,10,11,12};
Line Loop(3) = {2,8,-9,-7};

// DOMAIN SURFACE
Plane Surface(1) = {2};
// ABC SURFACE
Plane Surface(2) = {1};
Plane Surface(3) = {3};

// Define the rest of the surfaces
Physical Surface(0) = {1,2,3};
