// Gmsh project created on Thu Sep 22 08:50:48 2022
SetFactory("OpenCASCADE");
//+
lc = 2.0;
lc2 = 0.2;
Point(1) = {0, 0, 0, lc};
//+
Point(2) = {20, 0, 0, lc};
//+
Point(3) = {20, 20, 0, lc};
//+
Point(4) = {0, 20, 0, lc};
//+
Point(5) = {10, 10, 0, lc2};
//+
Line(1) = {4, 1};
//+
Line(2) = {1, 2};
//+
Line(3) = {2, 3};
//+
Line(4) = {3, 4};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Point{5} In Surface{1};
