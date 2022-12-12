// Gmsh project created on Mon Sep 05 15:26:03 2022
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
Plane Surface(1) = {1};
//+
Transfinite Curve {4, 3, 2, 1} = 10 Using Progression 1;
//+
Transfinite Surface {1};
//+
Transfinite Curve {4, 3, 2, 1} = 11 Using Progression 1;
//+
Transfinite Surface {1};
