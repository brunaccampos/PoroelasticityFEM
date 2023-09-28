// Gmsh project created on Wed Sep 27 11:25:57 2023
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 1, 0, 1.0};
//+
Point(3) = {0, 1, 0, 1.0};
//+
Point(4) = {1, 0, 0, 1.0};
//+
Point(5) = {0.22, 0, 0, 1.0};
//+
Point(6) = {0, 0.22, 0, 1.0};
//+
Circle(1) = {5, 1, 6};
//+
Line(2) = {4, 2};
//+
Line(3) = {5, 4};
//+
Line(4) = {2, 3};
//+
Line(5) = {3, 6};
//+
Curve Loop(1) = {5, -1, 3, 2, 4};
//+
Plane Surface(1) = {1};
//+
Transfinite Curve {1} = 21 Using Progression 1;
//+
Transfinite Curve {3, 2, 4, 5} = 11 Using Progression 1;
