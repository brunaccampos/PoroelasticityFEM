// Gmsh project created on Thu May 16 09:24:46 2024
SetFactory("OpenCASCADE");
// Points
Point(1) = {0, 0, 0, 1.0};
Point(2) = {37, 0, 0, 1.0};
Point(3) = {37, 7.3, 0, 1.0};
Point(4) = {27, 7.3, 0, 1.0};
Point(5) = {16, 14, 0, 1.0};
Point(6) = {12, 14, 0, 1.0};
Point(7) = {0, 14, 0, 1.0};
// Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 1};
// Curve
Transfinite Curve {1, 2, 3, 4, 5, 6, 7} = 10 Using Progression 1;
//
Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7};
Plane Surface(1) = {1};
