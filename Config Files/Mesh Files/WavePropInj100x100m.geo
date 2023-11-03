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

// New points 100m
Point(14) = {50, 0, 0, 1.0};
Point(15) = {75, 0, 0, 1.0};
Point(16) = {100, 0, 0, 1.0};
Point(17) = {100, 100, 0, 1.0};
Point(18) = {0, 100, 0, 1.0};
Point(19) = {0, 75, 0, 1.0};
Point(20) = {0, 50, 0, 1.0};
Point(21) = {75, 75, 0, 1.0};
Point(22) = {50, 50, 0, 1.0};

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
Line(18) = {2, 14};
Line(19) = {14, 15};
Line(20) = {15, 16};
Line(21) = {16, 17};
Line(22) = {17, 18};
Line(23) = {18, 19};
Line(24) = {19, 20};
Line(25) = {20, 4};
Line(26) = {19, 21};
Line(27) = {20, 22};
Line(28) = {21, 15};
Line(29) = {22, 14};
Line(30) = {17, 21};
Line(31) = {21, 22};
Line(32) = {22, 3};

// Delete auxiliar point
Delete {
  Point{13}; 
}

// Ratios
Nx1 = 31; Rx1 = 1.00;
Nx2 = 31; Rx2 = 1.00;
Nx3 = 21; Rx3 = 1.00;
Nx4 = 61; Rx4 = 0.95;
Nx5 = 51; Rx5 = 0.95;
Nx6 = 21; Rx6 = 1.00;
Nx7 = 21; Rx7 = 1.00;

// Curves
Transfinite Curve {6, 15, -3} = Nx1 Using Progression Rx1;
Transfinite Curve {5, 11, 9, 13, 4, 12, 10, 14, 21, 28, 29, 22, 26, 27} = Nx2 Using Progression Rx2;
Transfinite Curve {7, 16, 2} = Nx3 Using Progression Rx3;
Transfinite Curve {8, 17, 1} = Nx4 Using Progression Rx4;
Transfinite Curve {-18, 32, 25} = Nx5 Using Progression Rx5;
Transfinite Curve {19, 31, 24} = Nx6 Using Progression Rx6;
Transfinite Curve {20, 30, 23} = Nx7 Using Progression Rx7;

//+
Curve Loop(1) = {13, -17, -9, 8};
Curve Loop(2) = {14, 1, -10, 17};
Curve Loop(3) = {7, 9, -16, -11};
Curve Loop(4) = {10, 2, -12, 16};
Curve Loop(5) = {11, -15, 5, 6};
Curve Loop(6) = {12, 3, 4, 15};
Curve Loop(7) = {-5, -32, 29, -18};
Curve Loop(8) = {-29, -31, 28, -19};
Curve Loop(9) = {-28, -30, -21, -20};
Curve Loop(10) = {-25, 27, 32, -4};
Curve Loop(11) = {-24, 26, 31, -27};
Curve Loop(12) = {-23, -22, 30, -26};

//+
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};
Plane Surface(7) = {7};
Plane Surface(8) = {8};
Plane Surface(9) = {9};
Plane Surface(10) = {10};
Plane Surface(11) = {11};
Plane Surface(12) = {12};

//+
Transfinite Surface {1};
Transfinite Surface {2};
Transfinite Surface {3};
Transfinite Surface {4};
Transfinite Surface {5};
Transfinite Surface {6};
Transfinite Surface {7};
Transfinite Surface {8};
Transfinite Surface {9};
Transfinite Surface {10};
Transfinite Surface {11};
Transfinite Surface {12};