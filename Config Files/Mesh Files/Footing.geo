// Gmsh project created on Sun Sep 11 15:12:08 2022
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {8, 0, 0, 1.0};
//+
Point(3) = {8, -5, 0, 1.0};
//+
Point(4) = {0, -5, 0, 1.0};
//+
Point(5) = {1, 0, 0, 1.0};
//+
Line(1) = {1, 4};
//+
Line(2) = {4, 3};
//+
Line(3) = {3, 2};
//+
Line(4) = {2, 5};
//+
Line(5) = {5, 1};
//+
Curve Loop(1) = {1, 2, 3, 4, 5};
//+
Curve Loop(2) = {1, 2, 3, 4, 5};
//+
Plane Surface(1) = {2};
//+
MeshSize {5} = 0.1;
//+
MeshSize {1} = 0.1;
