function [Material, Mesh, BC, Control] = ColumnConsolidation1D_Dynamic()
% Column Consolidation 1D simulation
% Configuration File
% Based on Zienkiewicz (1982) model
%
% Assumptions/conventions:
% - stress is positive for tension
% - boundary condition for force is based on total stress
% - only solid acceleration is considered (undrained condition; no motions
% of the fluid relative to the solid skeleton can occur)
% - solid grains and fluid are incompressible

%% Material properties - Komijani (2019)

% elasticity modulus [Pa]
Material.E = 14.516e6;
% Poisson's ratio
Material.nu = 0.3;
% porous media permeability [m2/Pa s]
Material.kf = 1.0194e-6;
% dynamic viscosity [Pa s]
Material.mu = 1e-3;
% intrinsic permeability [m2]
Material.k = Material.kf * Material.mu;
% fluid bulk modulus [Pa]
Material.Kf = 2.1e9;
% solid bulk modulus [Pa]
Material.Ks = 1e20;
% material porosity
Material.n = 0.3;
% Biot's coefficient
Material.alpha = 1;
% fluid density [kg/m3]
Material.rho_f = 1000;
% solid density [kg/m3]
Material.rho_s = 2000;
% average density of the medium
Material.rho = Material.n*Material.rho_f + (1-Material.n)*Material.rho_s;
% 1/Q (related to storage coefficient)
Material.Q = (Material.alpha - Material.n)/Material.Ks + Material.n/Material.Kf;

%% Mesh parameters - uniform 1D mesh

% number of elements
ne = 10;

% column size [m]
L = 10;

% build mesh
Mesh = BuildMesh(ne, L);

%% Boundary conditions

% find top and bottom nodes for displacement field
BC.top_node_u = find(Mesh.coords_u == min(Mesh.coords_u));
BC.bottom_node_u = find(Mesh.coords_u == max(Mesh.coords_u));

% find top and bottom nodes for pressure field
BC.top_node_p = find(Mesh.coords_p == min(Mesh.coords_p));
BC.bottom_node_p = find(Mesh.coords_p == max(Mesh.coords_p));

% fixed nodes
%   displacement u=0 at the bottom
%   pressure p=0 at the top
%   no pressure gradient fp=0 at the bottom
BC.fixed_u = (BC.bottom_node_u);
BC.fixed_p = (BC.top_node_p);
BC.fixed_fp = (BC.bottom_node_p);

% free nodes
BC.free_u = setdiff(Mesh.dof_u, BC.fixed_u);
BC.free_p = setdiff(Mesh.dof_p, BC.fixed_p);
BC.free_fp = setdiff(Mesh.dof_p, BC.fixed_fp);

%% Dirichlet BCs
BC.fixed_u_value = zeros(length(BC.fixed_u),1);
BC.fixed_p_value = zeros(length(BC.fixed_p),1);

%% Neumann BCs
% applied traction [N]
BC.traction = 300;

BC.fext = zeros(Mesh.ndof_u, 1);
BC.fext(1,1) = BC.traction;

%% Quadrature order
Control.nq = 2;

%% Problem type
% 1 = steady state problem (no solid velocity, acceleration, and pressure
% change)
% 0 = transient problem (velocity and acceleration included)
Control.steady = 0;

%% Solution parameters
Control.dt = 5e-2;  % time step
Control.tend = 10;   % final simulation time
Control.tol = 1e-3; % tolerance for NR method
Control.max_it = 100; % maximum of iterations

% Control.plotu = round(length(Mesh.coords_u)/2);
% Control.plotp = round(length(Mesh.coords_p)/2);
Control.plotu = find(Mesh.coords_u == 6); % x = 6m
Control.plotp = find(Mesh.coords_p == 6); % x = 6m

%% Time discretization parameters
% Newmark method
Control.beta = 0.7;
Control.gamma = 0.7;
Control.theta = 0.7;

end