// Gmsh project created on Wed Sep 21 16:37:40 2022
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {120, 0, 0, 1.0};
//+
Recursive Delete {
  Point{2}; 
}
//+
Point(2) = {20, 0, 0, 1.0};
//+
Point(3) = {20, 20, 0, 1.0};
//+
Point(4) = {0, 20, 0, 1.0};
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
Transfinite Curve {1, 3, 4, 2} = 51 Using Progression 1;
//+
Transfinite Surface {1};
//+
Point(5) = {10, 10, 0, 1.0};
//+
Point(6) = {-3.5, 16.3, 0, 1.0};
//+
Recursive Delete {
  Point{6}; 
}
