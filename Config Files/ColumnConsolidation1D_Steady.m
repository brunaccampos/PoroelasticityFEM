function [Material, Mesh, BC, Control] = ColumnConsolidation1D_Steady()
% Column Consolidation 1D simulation

% Based on Korsawe (2006) model

% Assumptions/conventions:
% - stress is positive for tension
% - boundary condition for force is based on total stress
% - no acceleration terms for solid or fluid
% - solid velocity is neglected
% - fluid and solid grains are incompressible

%% Material properties - Korsawe (2006)

% elasticity modulus [Pa]
Material.E = 3e4;
% Poisson's ratio
Material.nu = 0.2;
% intrinsic permeability [m2]
Material.k = 1e-10;
% dynamic viscosity [Pa s]
Material.mu = 1e-3;
% porous media permeability [m2/Pa s]
Material.kf = Material.k/Material.mu;
% Biot's coefficient
Material.alpha = 1;
% 1/Q (related to storage coefficient)
Material.Q = 0; 

%% Mesh parameters - uniform 1D mesh

% number of elements
ne = 10;

% column size [m]
L = 1;

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
BC.traction = 1e3;

BC.fext = zeros(Mesh.ndof_u, 1);
BC.fext(1,1) = BC.traction;

%% Quadrature order
Control.nq = 2;

%% Problem type
% 1 = steady state problem
% 0 = transient problem
Control.steady = 1;

%% Solution parameters
Control.dt = 1;  % time step
Control.tend = 10;   % final simulation time
Control.plotu = round(length(Mesh.coords_u)/2);
Control.plotp = round(length(Mesh.coords_p)/2);

end