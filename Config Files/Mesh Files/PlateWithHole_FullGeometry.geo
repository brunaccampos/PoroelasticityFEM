// Gmsh project created on Mon Nov 14 14:27:22 2022
SetFactory("OpenCASCADE");

// Points
Point(1) = {0, 0, 0, 1.0};
Point(2) = {10, 0, 0, 1.0};
Point(3) = {10, 10, 0, 1.0};
Point(4) = {0, 10, 0, 1.0};
Point(5) = {5, 0, 0, 1.0};
Point(6) = {10, 5, 0, 1.0};
Point(7) = {5, 10, 0, 1.0};
Point(8) = {0, 5, 0, 1.0};
Point(9) = {5.5, 5, 0, 1.0};
Point(10) = {5.353553390593, 5.353553390593, 0, 1.0};
Point(11) = {5, 5.5, 0, 1.0};
Point(12) = {4.646446609406, 5.353553390593, 0, 1.0};
Point(13) = {4.5, 5, 0, 1.0};
Point(14) = {4.646446609406, 4.646446609406, 0, 1.0};
Point(15) = {5, 4.5, 0, 1.0};
Point(16) = {5.353553390593, 4.646446609406, 0, 1.0};

// Auxiliar point
Point(17) = {5, 5, 0, 1.0};

// Lines
Line(1) = {1, 5};
Line(2) = {5, 2};
Line(3) = {2, 6};
Line(4) = {6, 3};
Line(5) = {3, 7};
Line(6) = {7, 4};
Line(7) = {4, 8};
Line(8) = {8, 1};
Line(9) = {5, 15};
Line(10) = {2, 16};
Line(11) = {6, 9};
Line(12) = {3, 10};
Line(13) = {7, 11};
Line(14) = {4, 12};
Line(15) = {8, 13};
Line(16) = {1, 14};

// Circle
Circle(17) = {9, 17, 10};
Circle(18) = {10, 17, 11};
Circle(19) = {11, 17, 12};
Circle(20) = {12, 17, 13};
Circle(21) = {13, 17, 14};
Circle(22) = {14, 17, 15};
Circle(23) = {15, 17, 16};
Circle(24) = {16, 17, 9};

// Delete auxiliar point
Delete {
  Point{17};
}

// Ratios
Nx1 = 21; Rx1 = 1.00;
Nx2 = 41; Rx2 = 0.95;

// Curves
Transfinite Curve {1, 22} = Nx1 Using Progression Rx1;
Transfinite Curve {2, 23} = Nx1 Using Progression Rx1;
Transfinite Curve {3, 24} = Nx1 Using Progression Rx1;
Transfinite Curve {4, 17} = Nx1 Using Progression Rx1;
Transfinite Curve {5, 18} = Nx1 Using Progression Rx1;
Transfinite Curve {6, 19} = Nx1 Using Progression Rx1;
Transfinite Curve {7, 20} = Nx1 Using Progression Rx1;
Transfinite Curve {8, 21} = Nx1 Using Progression Rx1;
Transfinite Curve {9, 10, 11, 12, 13, 14, 15, 16} = Nx2 Using Progression Rx2;

// Curve loops
Curve Loop(1) = {9, 23, 10, 2};
Plane Surface(1) = {1};
Transfinite Surface {1};

Curve Loop(2) = {10, 24, 11, 3};
Plane Surface(2) = {2};
Transfinite Surface {2};

Curve Loop(3) = {11, 17, 12, 4};
Plane Surface(3) = {3};
Transfinite Surface {3};

Curve Loop(4) = {12, 18, 13, 5};
Plane Surface(4) = {4};
Transfinite Surface {4};

Curve Loop(5) = {13, 19, 14, 6};
Plane Surface(5) = {5};
Transfinite Surface {5};

Curve Loop(6) = {14, 20, 15, 7};
Plane Surface(6) = {6};
Transfinite Surface {6};

Curve Loop(7) = {15, 21, 16, 8};
Plane Surface(7) = {7};
Transfinite Surface {7};

Curve Loop(8) = {16, 22, 9, 1};
Plane Surface(8) = {8};
Transfinite Surface {8};