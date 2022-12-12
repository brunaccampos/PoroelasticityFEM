// Gmsh project created on Mon Dec 05 17:13:18 2022
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
Circle(5) = {7.5, 6.5, 0, 0.111125, 0, 2*Pi};
//+
Circle(6) = {7.5, 8.5, 0, 0.111125, 0, 2*Pi};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Curve Loop(2) = {5};
//+
Curve Loop(3) = {6};
//+
Plane Surface(1) = {1, 2, 3};