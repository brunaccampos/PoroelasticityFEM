// Gmsh project created on Thu Aug 18 15:24:19 2022
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {120, 0, 0, 1.0};
//+
Point(3) = {120, 100, 0, 1.0};
//+
Point(4) = {0, 100, 0, 1.0};
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
Recursive Delete {
  Surface{1}; 
}
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {120, 0, 0, 1.0};
//+
Point(3) = {120, 100, 0, 1.0};
//+
Point(4) = {0, 100, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Circle(5) = {60, 40, 0, 0.1, 0, 2*Pi};
//+
Circle(6) = {60, 35, 0, 0.1, 0, 2*Pi};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Curve Loop(2) = {5};
//+
Curve Loop(3) = {6};
//+
Plane Surface(1) = {1, 2, 3};
//+
Transfinite Surface {1};
//+
Transfinite Curve {4, 2} = 100 Using Progression 1;
//+
Transfinite Curve {1, 3} = 120 Using Progression 1;
//+
Transfinite Surface {1};
//+
Transfinite Curve {1, 3} = 60 Using Progression 1;
//+
Transfinite Curve {4, 2} = 50 Using Progression 1;
//+
Point(7) = {0, 25, 0, 1.0};
//+
Point(8) = {120, 25, 0, 1.0};
//+
Point(9) = {120, 50, 0, 1.0};
//+
Point(10) = {0, 50, 0, 1.0};
//+
Line(7) = {10, 9};
//+
Line(8) = {7, 8};
//+
Curve Loop(4) = {4, 1, 2, 3};
//+
Recursive Delete {
  Curve{4}; Curve{2}; 
}
//+
Recursive Delete {
  Curve{7}; Surface{1}; 
}
//+
Recursive Delete {
  Curve{8}; 
}
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {120, 0, 0, 1.0};
//+
Point(3) = {120, 100, 0, 1.0};
//+
Point(4) = {0, 100, 0, 1.0};
//+
Point(5) = {0, 25, 0, 1.0};
//+
Point(6) = {0, 50, 0, 1.0};
//+
Point(7) = {120, 50, 0, 1.0};
//+
Point(8) = {120, 25, 0, 1.0};
//+
Line(1) = {4, 6};
//+
Line(2) = {6, 7};
//+
Line(3) = {7, 3};
//+
Line(4) = {3, 4};
//+
Line(5) = {6, 5};
//+
Line(6) = {5, 8};
//+
Line(7) = {8, 7};
//+
Line(8) = {5, 1};
//+
Line(9) = {1, 2};
//+
Line(10) = {2, 8};
//+
Circle(11) = {60, 35, 0, 0.1, 0, 2*Pi};
//+
Circle(12) = {60, 40, 0, 0.1, 0, 2*Pi};
//+
Curve Loop(5) = {1, 2, 3, 4};
//+
Plane Surface(1) = {5};
//+
Curve Loop(6) = {5, 6, 7, -2};
//+
Curve Loop(7) = {12};
//+
Curve Loop(8) = {11};
//+
Plane Surface(2) = {6, 7, 8};
//+
Curve Loop(9) = {9, 10, -6, 8};
//+
Plane Surface(3) = {9};
//+
Transfinite Surface {1};
//+
Transfinite Surface {3};
//+
Transfinite Surface {2};
//+
Transfinite Curve {8, 10, 3, 1} = 10 Using Progression 1;
//+
Transfinite Curve {4, 2, 9, 6} = 20 Using Progression 1;
//+
Transfinite Curve {5, 7} = 10 Using Progression 1;
//+
Curve Loop(10) = {1, 2, 3, 4};
//+
Plane Surface(4) = {10};
//+
Curve Loop(11) = {5, 6, 7, -2};
//+
Curve Loop(12) = {12};
//+
Curve Loop(13) = {11};
//+
Plane Surface(5) = {11, 12, 13};
//+
Curve Loop(14) = {9, 10, -6, 8};
//+
Plane Surface(6) = {14};
//+
Transfinite Surface {4};
//+
Transfinite Surface {5};
//+
Transfinite Surface {6};
//+
Transfinite Curve {8, 10, 5, 7} = 10 Using Progression 1;
//+
Transfinite Curve {1, 3} = 20 Using Progression 1;
//+
Transfinite Curve {4, 2, 6, 9} = 25 Using Progression 1;
