// Gmsh project created on Thu Oct 26 14:59:09 2023
SetFactory("OpenCASCADE");
// Points
Point(1) = {0.22, 0, 0, 1.0};
Point(2) = {10, 0, 0, 1.0};
Point(3) = {10, 10, 0, 1.0};
Point(4) = {0, 10, 0, 1.0};
Point(5) = {0, 0.22, 0, 1.0};
Point(6) = {2.5, 0, 0, 1.0};
Point(7) = {2.5, 2.5, 0, 1.0};
Point(8) = {0, 2.5, 0, 1.0};
Point(9) = {5, 0, 0, 1.0};
Point(10) = {5, 5, 0, 1.0};
Point(11) = {0, 5, 0, 1.0};
Point(12) = {0.155563491861, 0.155563491861, 0, 1.0};

// Auxiliar point
Point(13) = {0, 0, 0, 1.0};

// Lines
Line(1) = {8, 5};
Line(2) = {8, 11};
Line(3) = {11, 4};
Line(4) = {4, 3};
Line(5) = {3, 2};
Line(6) = {2, 9};
Line(7) = {9, 6};
Line(8) = {6, 1};
Line(9) = {6, 7};
Line(10) = {7, 8};
Line(11) = {9, 10};
Line(12) = {10, 11};
Circle(13) = {1, 13, 12};
Circle(14) = {12, 13, 5};
Line(15) = {3, 10};
Line(16) = {10, 7};
Line(17) = {7, 12};

// Delete auxiliar point
Delete {
  Point{13}; 
}

// Ratios
Nx1 = 21; Rx1 = 1.00;
Nx2 = 21; Rx2 = 1.00;
Nx3 = 21; Rx3 = 1.00;
Nx4 = 51; Rx4 = 0.95;

// Curves
Transfinite Curve {6, 15, 3} = Nx1 Using Progression Rx1;
Transfinite Curve {5, 11, 9, 13, 4, 12, 10, 14} = Nx2 Using Progression Rx2;
Transfinite Curve {7, 16, 2} = Nx3 Using Progression Rx3;
Transfinite Curve {8, 17, 1} = Nx4 Using Progression Rx4;

//+
Curve Loop(1) = {13, -17, -9, 8};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {14, 1, -10, 17};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {7, 9, -16, -11};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {10, 2, -12, 16};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {11, -15, 5, 6};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {12, 3, 4, 15};
//+
Plane Surface(6) = {6};
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Transfinite Surface {3};
//+
Transfinite Surface {4};
//+
Transfinite Surface {5};
//+
Transfinite Surface {6};
