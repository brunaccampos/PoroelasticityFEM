// Gmsh project created on Fri Aug 30 15:38:50 2024
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {4000, 0, 0, 1.0};
//+
Point(3) = {4000, 800, 0, 1.0};
//+
Point(4) = {4000, 1200, 0, 1.0};
//+
Point(5) = {4000, 2000, 0, 1.0};
//+
Point(6) = {0, 2000, 0, 1.0};
//+
Point(7) = {0, 1200, 0, 1.0};
//+
Point(8) = {0, 800, 0, 1.0};
//+
Point(9) = {1000, 0, 0, 1.0};
//+
Point(10) = {1000, 800, 0, 1.0};
//+ 
Point(11) = {1000, 1200, 0, 1.0};
//+
Point(12) = {1000, 2000, 0, 1.0};
//+ 
Point(13) = {3000, 0, 0, 1.0};
//+
Point(14) = {3000, 800, 0, 1.0};
//+
Point(15) = {3000, 1200, 0, 1.0};
//+
Point(16) = {3000, 2000, 0, 1.0};
//+
Line(1) = {1, 9};
//+
Line(2) = {9, 13};
//+
Line(3) = {13, 2};
//+
Line(4) = {2, 3};
//+
Line(5) = {3, 4};
//+
Line(6) = {4, 5};
//+
Line(7) = {5, 16};
//+
Line(8) = {16, 12};
//+ 
Line(9) = {12, 6};
//+
Line(10) = {6, 7};
//+
Line(11) = {7, 8};
//+
Line(12) = {8, 1};
//+
Line(13) = {7, 11};
//+
Line(14) = {11, 15};
//+
Line(15) = {15, 4};
//+
Line(16) = {8, 10};
//+
Line(17) = {10, 14};
//+
Line(18) = {14, 3};
//+
Line(19) = {12, 11};
//+
Line(20) = {11, 10};
//+
Line(21) = {10, 9};
//+
Line(22) = {16, 15};
//+
Line(23) = {15, 14};
//+
Line(24) = {14, 13};

// Ratios
N1 = 151; R1 = 1.00;
N2 = 31; R2 = 0.98;
N3 = 21; R3 = 1.00;
N4 = 41; R4 = 0.98;

// Curves
Transfinite Curve {2, 17, 14, -8} = N1 Using Progression R1;
Transfinite Curve {1, 16, 13, -9, -3, -18, -15, 7} = N2 Using Progression R2;
Transfinite Curve {11, 20, 23, -5} = N3 Using Progression R3;
Transfinite Curve {-12, -21, -24, 4, 10, 19, 22, -6} = N4 Using Progression R4;

// 
Curve Loop(1) = {1, -21, -16, 12};
Plane Surface(1) = {1};
// 
Curve Loop(2) = {2, -24, -17, 21};
Plane Surface(2) = {2};
//
Curve Loop(3) = {3, 4, -18, 24};
Plane Surface(3) = {3};
//
Curve Loop(4) = {16, -20, -13, 11};
Plane Surface(4) = {4};
//
Curve Loop(5) = {17, -23, -14, 20};
Plane Surface(5) = {5};
//
Curve Loop(6) = {18, 5, -15, 23};
Plane Surface(6) = {6};
//
Curve Loop(7) = {13, -19, 9, 10};
Plane Surface(7) = {7};
//
Curve Loop(8) = {14, -22, 8, 19};
Plane Surface(8) = {8};
//
Curve Loop(9) = {15, 6, 7, 22};
Plane Surface(9) = {9};

//
Transfinite Surface {1};
Transfinite Surface {2};
Transfinite Surface {3};
Transfinite Surface {4};
Transfinite Surface {5};
Transfinite Surface {6};
Transfinite Surface {7};
Transfinite Surface {8};
Transfinite Surface {9};
