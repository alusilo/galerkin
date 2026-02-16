//=============================
// Point source geometry
//=============================

// Fault dimension
f_w = 3000.; // width (x - coordinate)

// Nucleation dimension
n_w = 300.;  // width

// Length of the simulation zone
dx_sz = 1000.;
dy_sz = 1500.;

// Length of the CPML region
dx_cpml = 300.;
dy_cpml = 300.;

// Characteristic length for
// * simulation area
len = 100.;
// * CPML region
len_cpml = 100.;

// To avoid negative coordinates
min_x = dx_sz+dx_cpml+(f_w/2.0);
min_y = dy_sz+dy_cpml;

// define the 2 points of the nucleation patch
Point(1) = {min_x-(n_w/2.0), min_y, 0.0, len};
Point(2) = {min_x+(n_w/2.0), min_y, 0.0, len};
// define the 2 points of the right patch
Point(3) = {min_x-(n_w/2.0)+(f_w/4.0), min_y, 0.0, len};
Point(4) = {min_x+(n_w/2.0)+(f_w/4.0), min_y, 0.0, len};
// define the 2 points of the left patch
Point(5) = {min_x-(n_w/2.0)-(f_w/4.0), min_y, 0.0, len};
Point(6) = {min_x+(n_w/2.0)-(f_w/4.0), min_y, 0.0, len};
// define the limits of the fault
Point(7) = {min_x-(f_w/2.0), min_y, 0.0, len};
Point(8) = {min_x+(f_w/2.0), min_y, 0.0, len};
// define the 4 points of the lower simulation zone
Point(9) = {min_x-dx_sz-(f_w/2.0), min_y-dy_sz, 0.0, len_cpml};
Point(10) = {min_x+dx_sz+(f_w/2.0), min_y-dy_sz, 0.0, len_cpml};
Point(11) = {min_x+dx_sz+(f_w/2.0), min_y, 0.0, len_cpml};
Point(12) = {min_x-dx_sz-(f_w/2.0), min_y, 0.0, len_cpml};
// define the 2 points of the upper simulation zone
Point(13) = {min_x+dx_sz+(f_w/2.0), min_y+dy_sz, 0.0, len_cpml};
Point(14) = {min_x-dx_sz-(f_w/2.0), min_y+dy_sz, 0.0, len_cpml};
// define the 4 points of the CPML region
Point(15) = {min_x-dx_sz-(f_w/2.0)-dx_cpml, min_y-dy_sz-dy_cpml, 0.0, len_cpml};
Point(16) = {min_x+dx_sz+(f_w/2.0)+dx_cpml, min_y-dy_sz-dy_cpml, 0.0, len_cpml};
Point(17) = {min_x+dx_sz+(f_w/2.0)+dx_cpml, min_y+dy_sz+dy_cpml, 0.0, len_cpml};
Point(18) = {min_x-dx_sz-(f_w/2.0)-dx_cpml, min_y+dy_sz+dy_cpml, 0.0, len_cpml};
// Define the physical points
Physical Point(0) = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18};

// define the edge of the nucleation patch
Line(1) = {1,2};
// define the edge of the right patch
Line(2) = {3,4};
// define the edge of the left patch
Line(3) = {5,6};
// define the rest of the fault surface
Line(4) = {7,5};
Line(5) = {6,1};
Line(6) = {2,3};
Line(7) = {4,8};
// define the edges of the lower simulation zone
Line(8) = {9,10};
Line(9) = {10,11};
Line(10) = {11,8};
Line(11) = {7,12};
Line(12) = {12,9};
// define the edges of the upper simulation zone 
Line(13) = {11,13};
Line(14) = {13,14};
Line(15) = {14,12};
// define the 4 edges of the CPML region
Line(16) = {15,16};
Line(17) = {16,17};
Line(18) = {17,18};
Line(19) = {18,15};

// Define the physical lines
Physical Line(0) = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};

// define the surface of the lower simulation zone
Line Loop(1) = {8,9,10,-7,-2,-6,-1,-5,-3,-4,11,12};
Plane Surface(1) = {1};
// define the surface of the upper simulation zone
Line Loop(2) = {-11,4,3,5,1,6,2,7,-10,13,14,15};
Plane Surface(2) = {2};
// define the surface of the CPML region
Line Loop(3) = {16,17,18,19};
Line Loop(4) = {8,9,13,14,15,12};
Plane Surface(3) = {3,4};
// Define the rest of the surfaces
Physical Surface(0) = {1,2,3};
