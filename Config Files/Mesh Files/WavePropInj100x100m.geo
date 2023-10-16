// Gmsh project created on Wed Sep 27 11:25:57 2023
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 10.0};
//+
Point(2) = {100, 100, 0, 10.0};
//+
Point(3) = {0, 100, 0, 10.0};
//+
Point(4) = {100, 0, 0, 10.0};
//+
Point(5) = {0.22, 0, 0, 0.01};
//+
Point(6) = {0, 0.22, 0, 0.01};
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
Recursive Delete {
  Point{1}; 
}
