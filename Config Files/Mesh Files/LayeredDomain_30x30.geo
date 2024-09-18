// Gmsh project created on Fri Aug 30 15:38:50 2024
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {30, 0, 0, 1.0};
//+
Point(3) = {30, 10, 0, 1.0};
//+
Point(4) = {30, 30, 0, 1.0};
//+
Point(5) = {0, 30, 0, 1.0};
//+
Point(6) = {0, 10, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 1};
//+
Line(7) = {6, 3};

// Ratios
N1 = 101; R1 = 1.00;
N2 = 101; R2 = 1.00;

// Curves
Transfinite Curve {1, 7, 4} = N1 Using Progression R1;
Transfinite Curve {-6, 2, 5, -3} = N2 Using Progression R2;

// 
Curve Loop(1) = {1, 2, -7, 6};
Plane Surface(1) = {1};
// 
Curve Loop(2) = {7, 3, 4, 5};
Plane Surface(2) = {2};

//
Transfinite Surface {1};
Transfinite Surface {2};
