// Gmsh project created on Fri Aug 30 15:38:50 2024
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {30, 0, 0, 1.0};
//+
Point(3) = {30, 10, 0, 1.0};
//+
Point(4) = {30, 20, 0, 1.0};
//+
Point(5) = {30, 30, 0, 1.0};
//+
Point(6) = {0, 30, 0, 1.0};
//+
Point(7) = {0, 20, 0, 1.0};
//+
Point(8) = {0, 10, 0, 1.0};
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
Line(6) = {6, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 1};
//+
Line(9) = {7, 4};
//+
Line(10) = {8, 3};

// Ratios
N1 = 141; R1 = 1.00;
N2 = 41; R2 = 1.00;
N3 = 41; R3 = 1.00;

// Curves
Transfinite Curve {1, 5, 9, 10} = N1 Using Progression R1;
Transfinite Curve {2, 4, 6, 8} = N2 Using Progression R2;
Transfinite Curve {3, 7} = N3 Using Progression R3;

// 
Curve Loop(1) = {6, 9, 4, 5};
Plane Surface(1) = {1};
// 
Curve Loop(2) = {7, 10, 3, -9};
Plane Surface(2) = {2};
//
Curve Loop(3) = {8, 1, 2, -10};
Plane Surface(3) = {3};

//
Transfinite Surface {1};
Transfinite Surface {2};
Transfinite Surface {3};
