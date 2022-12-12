// Gmsh project created on Wed Jun 22 12:59:51 2022
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0.1, 0, 0, 1.0};
//+
Point(3) = {0.1, 10, 0, 1.0};
//+
Point(4) = {0, 10, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Plane Surface(1) = {1};
//+
Transfinite Curve {3, 2, 1, 4} = 10 Using Progression 1;
//+
Transfinite Surface {1};
//+
Recombine Surface {1};
