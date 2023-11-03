// Gmsh project created on Thu Oct 26 14:59:09 2023
SetFactory("OpenCASCADE");
// Points
Point(1) = {0.22, 0, 0, 1.0};
Point(2) = {0, 0.22, 0, 1.0};
Point(3) = {0.155563491861, 0.155563491861, 0, 1.0};
Point(4) = {100, 0, 0, 1.0};
Point(5) = {100, 100, 0, 1.0};
Point(6) = {0, 100, 0, 1.0};

// Auxiliar point
Point(7) = {0, 0, 0, 1.0};

// Lines
Line(1) = {1, 4};
Line(2) = {4, 5};
Line(3) = {5, 6};
Line(4) = {2, 6};
Line(5) = {3, 5};
Circle(6) = {1, 7, 3};
Circle(7) = {3, 7, 2};

// Delete auxiliar point
Delete {
  Point{7}; 
}

// Ratios
Nx1 = 251; Rx1 = 0.97;
Nx2 = 41; Rx2 = 1.00;

// Curves
Transfinite Curve {-1, -5, -4} = Nx1 Using Progression Rx1;
Transfinite Curve {2, 6, 3, 7} = Nx2 Using Progression Rx2;

//+
Curve Loop(1) = {6, 5, -2, -1};
Curve Loop(2) = {7, 4, -3, -5};

//+
Plane Surface(1) = {1};
Plane Surface(2) = {2};

//+
Transfinite Surface {1};
Transfinite Surface {2};