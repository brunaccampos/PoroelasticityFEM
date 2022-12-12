function [Material, MeshU, MeshP, MeshN, BC, Control] = Column1D_Dynamic_Komijani(config_dir, progress_on)
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

%% Poroelasticity model
% 1 - Biot theory
% 0 - Spanos theory (additional porosity equation)
Control.Biotmodel = 1;

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
Material.Qinv = (Material.alpha - Material.n)/Material.Ks + Material.n/Material.Kf;

%% Mesh parameters
if progress_on
    disp([num2str(toc),': Building Mesh...']);
end

% mesh type
% 'Manual': 1D mesh
% 'Gmsh': 2D mesh, input file from GMSH
MeshType = 'Manual';

switch MeshType
    case 'Manual'
        % number of space dimensions
        nsd = 1;
        % number of elements
        ne = 10;
        % column size [m]
        L = 10;
        %%%% solid displacement field
        typeU = 'L3';
        MeshU = Build1DMesh(nsd, ne, L, typeU);
        %%%% fluid pressure field
        typeP = 'L2';
        MeshP = Build1DMesh(nsd, ne, L, typeP);
        %%%% porosity field
        if ~Control.Biotmodel
            typeN = 'L2';
            MeshN = Build1DMesh(nsd, ne, L, typeN);
        else
            MeshN = [];
        end
    case 'Gmsh'
        % Version 2 ASCII
        % number of space dimensions
        nsd = 2;
        %%%% displacement field
        fieldU = 'u';
        meshFileNameU = 'Column2DQ9.msh';
        MeshU = BuildMesh_GMSH(meshFileNameU, fieldU, nsd, config_dir, progress_on);
        %%%% pressure field
        fieldP = 'p';
        meshFileNameP = 'Column2DQ4.msh';
        MeshP = BuildMesh_GMSH(meshFileNameP, fieldP, nsd, config_dir, progress_on);
        %%%% porosity field
        if ~Control.Biotmodel
            fieldN = 'n';
            meshFileNameN = 'Column2DQ4.msh';
            MeshN = BuildMesh_GMSH(meshFileNameN, fieldN, nsd, config_dir, progress_on);
        else
            MeshN = [];
        end
end

%% Boundary conditions

% find top and bottom nodes for displacement field
BC.top_node_u = find(MeshU.coords == min(MeshU.coords));
BC.bottom_node_u = find(MeshU.coords == max(MeshU.coords));

% find top and bottom nodes for pressure field
BC.top_node_p = find(MeshP.coords == min(MeshP.coords));
BC.bottom_node_p = find(MeshP.coords == max(MeshP.coords));

% fixed nodes
%   displacement u=0 at the bottom
%   pressure p=0 at the top
%   no pressure gradient fp=0 at the bottom
BC.fixed_u = (BC.bottom_node_u);
BC.fixed_p = (BC.top_node_p);
BC.fixed_fp = (BC.bottom_node_p);

% free nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);
BC.free_fp = setdiff(MeshP.DOF, BC.fixed_fp);

%% Dirichlet BCs
BC.fixed_u_value = zeros(length(BC.fixed_u),1);
BC.fixed_p_value = zeros(length(BC.fixed_p),1);

%% Neumann BCs
% point traction [N]
BC.pointLoadValue = 3000;
BC.pointLoadNodes = BC.top_node_u;
BC.pointLoad = zeros(MeshU.nDOF,1);
BC.pointLoad(BC.pointLoadNodes) = BC.pointLoadValue;

% distributed traction [N/m]
BC.tractionNodes = [];

% body force
BC.b = @(x)[];  

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

Control.plotu = find(MeshU.coords == 5); % x = 6m
Control.plotp = find(MeshP.coords == 5); % x = 6m

%% Time discretization parameters
% Newmark method
Control.beta = 0.7;
Control.gamma = 0.7;
Control.theta = 0.7;

end