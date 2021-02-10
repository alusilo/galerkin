//=============================
// Point source geometry
//=============================

x_min = 0.;
x_max = 2.; // max width
y_min = 0.;
y_max = 2.; // max heigth

// Characteristic length for
// * simulation area
vs = 1.0;
vp = 2.0;
lambda = 1.4142135624;
freq = vp/lambda;
len = vs/(3*freq);

Point(1) = {x_min, y_min, 0.0, len};
Point(2) = {x_min, y_max, 0.0, len};
Point(3) = {x_max, y_max, 0.0, len};
Point(4) = {x_max, y_min, 0.0, len};

// Define the physical points
Physical Point(0) = {1,2,3,4};

// define the edges
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

// Define the physical lines
Physical Line(0) = {1,2,3,4};

// define the surface
Line Loop(1) = {1,2,3,4};

Plane Surface(1) = {1};

// Define the rest of the surfaces
Physical Surface(0) = {1};
