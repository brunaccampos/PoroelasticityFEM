function [Material, MeshU, MeshP, MeshN, BC, Control] = Column2D_Dynamic_Komijani(config_dir, progress_on)
% Column Consolidation 2D simulation
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
MeshType = 'Gmsh';

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

%% Dirichlet BCs
% displacement u=0 at bottom (y), left (x), and right (x)
BC.fixed_u = [MeshU.left_dofx; MeshU.right_dofx; MeshU.bottom_dofy];
% pressure p=0 at top
BC.fixed_p = [MeshP.top_dof];

% fixed DOF values
BC.fixed_u_value = zeros(length(BC.fixed_u),1);
BC.fixed_p_value = zeros(length(BC.fixed_p),1);

% free nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs
% point loads
BC.pointLoad = [];

% impervious at bottom, left, and right
BC.fixed_fp = [MeshP.left_dof; MeshP.right_dof; MeshP.bottom_dof];
% free nodes
BC.free_fp = setdiff(MeshP.DOF, BC.fixed_fp);

% prescribed traction
BC.traction = 3000; % magnitute [N/m^2]
BC.tractionNodes = MeshU.top_nodes;
Force = BC.traction * max(MeshU.coords(:,1))/((length(MeshU.top_nodes) - 1)/2);
BC.tractionForce = zeros(length(BC.tractionNodes),2);

% Q9 elements for displacement field
for n = 1:length(BC.tractionForce)
    if any(BC.tractionNodes(n) == MeshU.conn(:,1:4),'all') % then node is a corner node
        BC.tractionForce(n,:) = [0, Force/3];
    else % then node is a midside node
        BC.tractionForce(n,:) = [0, Force*2/3];
    end
end

% find the nodes in the top left and right corners
lefttopnode = find(MeshU.coords(BC.tractionNodes,1) == min(MeshU.coords(:,1)));
righttopnode  = find(MeshU.coords(BC.tractionNodes,1) == max(MeshU.coords(:,1)));

BC.tractionForce(lefttopnode,2) = BC.tractionForce(lefttopnode,2)/2;
BC.tractionForce(righttopnode,2) = BC.tractionForce(righttopnode,2)/2;

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
Control.dt = 1e-2;  % time step
Control.tend = 10;   % final simulation time
Control.tol = 1e-3; % tolerance for NR method
Control.max_it = 100; % maximum of iterations

Control.plotu = 184; % dof y of node 92 (x = 0.05m, y = 5m)
Control.plotp = 52; % dof of node 52 (x = 0.05m, y = 5m)

%% Time discretization parameters
% Newmark method
Control.beta = 0.7;
Control.gamma = 0.7;
Control.theta = 0.7;

end