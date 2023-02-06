// Gmsh project created on Mon Feb 06 13:01:46 2023
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {15, 0, 0, 1.0};
//+
Point(3) = {15, 7, 0, 1.0};
//+
Point(4) = {0, 7, 0, 1.0};
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
Curve Loop(2) = {3, 4, 1, 2};
//+
Plane Surface(2) = {2};
//+
Transfinite Surface {2};
//+
Transfinite Curve {1, 3} = 31 Using Progression 1;
//+
Transfinite Curve {4, 2} = 15 Using Progression 1;
