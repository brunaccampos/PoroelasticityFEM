// Gmsh project created on Thu Oct 20 15:38:22 2022
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {10, 0, 0, 1.0};
//+
Point(3) = {10, 10, 0, 1.0};
//+
Point(4) = {0, 10, 0, 1.0};
//+
Point(5) = {5, 10, 0, 1.0};
//+
Line(1) = {4, 1};
//+
Line(2) = {1, 2};
//+
Line(3) = {2, 3};
//+
Line(4) = {3, 5};
//+
Line(5) = {5, 4};
//+
Curve Loop(1) = {1, 2, 3, 5, 4};
//+
Curve Loop(2) = {1, 2, 3, 5, 4};
//+
Plane Surface(1) = {2};
//+
MeshSize {5} = 0.1;

//+
Transfinite Surface {1};
