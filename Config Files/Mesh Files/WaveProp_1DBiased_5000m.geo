// Gmsh project created on Fri Oct 27 13:55:24 2023
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {2500, 0, 0, 1.0};
//+
Point(3) = {5000, 0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Transfinite Curve {1} = 251 Using Progression 0.95;
//+
Transfinite Curve {-2} = 251 Using Progression 0.95;
