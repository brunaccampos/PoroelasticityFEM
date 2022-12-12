// Gmsh project created on Fri Sep 09 13:41:01 2022
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {120, 0, 0, 1.0};
//+
Point(3) = {120, 25, 0, 1.0};
//+
Point(4) = {0, 25, 0, 1.0};
//+
Line(1) = {4, 1};
//+
Line(2) = {1, 2};
//+
Line(3) = {2, 3};
//+
Line(4) = {3, 4};
//+
Circle(5) = {60, 10, 0.1, 0.5, 0, 2*Pi};
//+
Recursive Delete {
  Curve{5}; 
}
//+
Circle(5) = {60, 15, 0, 0.1, 0, 2*Pi};
//+
Circle(6) = {60, 10, 0, 0.1, 0, 2*Pi};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Curve Loop(2) = {5};
//+
Curve Loop(3) = {6};
//+
Plane Surface(1) = {1, 2, 3};
