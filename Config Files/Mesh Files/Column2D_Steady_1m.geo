// Gmsh project created on Wed Aug 17 15:18:46 2022
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Curve Loop(2) = {4, 1, 2, 3};
//+
Plane Surface(1) = {2};
//+
Transfinite Surface {1};
//+
Transfinite Curve {4, 1, 2, 3} = 10 Using Progression 1;
