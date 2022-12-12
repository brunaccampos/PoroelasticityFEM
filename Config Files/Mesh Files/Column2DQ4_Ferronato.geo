// Gmsh project created on Thu Nov 10 09:31:07 2022
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0.1, 0, 0, 1.0};
//+
Point(3) = {0.1, 15, 0, 1.0};
//+
Point(4) = {0, 15, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Transfinite Curve {4, 2} = 31 Using Progression 1;
//+
Transfinite Curve {1, 3} = 4 Using Progression 1;
//+
Transfinite Surface {1};
