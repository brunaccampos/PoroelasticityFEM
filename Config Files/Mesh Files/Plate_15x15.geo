// Gmsh project created on Fri May 26 15:11:27 2023
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {15, 0, 0, 1.0};
//+
Point(3) = {15, 15, 0, 1.0};
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
Transfinite Surface {1};
//+
Transfinite Surface {1} = {1, 2, 3, 4};
