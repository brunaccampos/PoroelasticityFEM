function [Material, MeshU, MeshP, MeshN, BC, Control] = Footing2D_Diebels(config_dir, progress_on)
% 2D simulation of footing problem
% Configuration File
% Based on Korsawe (2006) model
% ------------------------------------------------------------------------
% Assumptions/conventions:
% - stress is positive for tension
% - boundary condition for force is based on total stress
% - no acceleration terms for solid or fluid
% - solid velocity is neglected
% - fluid and solid grains are incompressible
% - porosity is constant in space and varies over time
% ------------------------------------------------------------------------

%% Poroelasticity model
% 1 - Biot theory
% 0 - Spanos theory (additional porosity equation)
Control.Biotmodel = 0;

%% Material properties - Diebels (1996)
% Lame constant [GPa]
lambda = 5.583e-3;
% shear modulus [GPa]
G = 8.375e-3;
% elasticity modulus [GPa]
Material.E = G*(3*lambda+2*G)/(lambda+G);
% Poisson's ratio
Material.nu = 0.2;
% porous media permeability [m2/GPa s]
Material.kf = 1e3;
% dynamic viscosity [GPa s]
Material.mu = 1e-12;
% intrinsic permeability [m2]
Material.k = Material.kf * Material.mu;
% material porosity
Material.n = 0.33;
% Biot's coefficient
Material.alpha = 1;
% 1/Q (related to storage coefficient)
Material.Minv = 0;

% lumped mass matrix - 0: false, 1: true
Material.lumpedMass = 0;

% constititive law - 'PlaneStress' or 'PlaneStrain'
% Note: use 'PlaneStrain' for 1D or 2D poroelasticity
Material.constLaw = 'PlaneStrain';

%% Spanos material parameters
% porosity effective pressure coefficient (Spanos, 1989)
% n = 0; % lower limit
n = 1; % return to Biot
% n = Material.Ks/Material.Kf; % upper limit
 
% porosity equation coefficients
Material.deltaF = 0;
Material.deltaS = Material.alpha - Material.n;

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
        L = 1;
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
        meshFileNameU = 'Mesh Files\Footing_DiebelsQ9.msh';
        MeshU = BuildMesh_GMSH(meshFileNameU, fieldU, nsd, config_dir, progress_on);
        %%%% pressure field
        fieldP = 'p';
        meshFileNameP = 'Mesh Files\Footing_DiebelsQ4.msh';
        MeshP = BuildMesh_GMSH(meshFileNameP, fieldP, nsd, config_dir, progress_on);
        %%%% porosity field
        if ~Control.Biotmodel
            fieldN = 'n';
            meshFileNameN = 'Mesh Files\Footing_DiebelsQ4.msh';
            MeshN = BuildMesh_GMSH(meshFileNameN, fieldN, nsd, config_dir, progress_on);
        else
            MeshN = [];
        end
end

%% Initial conditions
% displacement
BC.initU = [];

% pressure
BC.initP = [];

%% Load strip 
% select nodes in the interval [(0,10) (5,10)] where the load is applied
coord1 = 5;
coord2 = 10;
nu = 1;
for i = 1:length(MeshU.top_nodes)
   if MeshU.coords(MeshU.top_nodes(i),1) >= coord1 &&  MeshU.coords(MeshU.top_nodes(i),1) <= coord2
      nodesU(nu,1) = MeshU.top_nodes(i);
      nu = nu+1;
   end
end

np = 1;
for i = 1:length(MeshP.top_nodes)
   if MeshP.coords(MeshP.top_nodes(i),1) >= coord1 &&  MeshP.coords(MeshP.top_nodes(i),1) <= coord2
      nodesP(np,1) = MeshP.top_nodes(i);
      np = np+1;
   end
end

%% Dirichlet BCs - solid
% displacement u=0 at bottom (y), left (x), and right (x)
BC.fixed_u = [MeshU.left_dofx; MeshU.right_dofx; MeshU.bottom_dofy];
% fixed DOF values
BC.fixed_u_value = zeros(length(BC.fixed_u),1);
% free nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid
% pressure p=0 at top
BC.fixed_p = setdiff(MeshP.top_nodes, nodesP);
% fixed DOF values
BC.fixed_p_value = zeros(length(BC.fixed_p),1);
% free nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% traction interpolation (needed for traction applied in wells); 1 - true, 0 - false
BC.tractionInterp = 0;
% prescribed traction [GN/m2]
BC.traction = -15e-6;
BC.tractionNodes = nodesU;
Force = BC.traction * 1/((length(BC.tractionNodes) - 1)/2);
BC.tractionForce = zeros(length(BC.tractionNodes),2);

% Q9 elements for displacement field
for n = 1:length(BC.tractionForce)
    if any(BC.tractionNodes(n) == MeshU.conn(:,1:4),'all') % then node is a corner node
        BC.tractionForce(n,:) = [0, Force/3];
    else % then node is a midside node
        BC.tractionForce(n,:) = [0, Force*2/3];
    end
end

% find the nodes in the left and right corners of the strip where load is
% applied
lefttopnode = find(MeshU.coords(BC.tractionNodes,1) == 5);
righttopnode  = find(MeshU.coords(BC.tractionNodes,1) == max(MeshU.coords(:,1)));

BC.tractionForce(lefttopnode,2) = BC.tractionForce(lefttopnode,2)/2;
BC.tractionForce(righttopnode,2) = BC.tractionForce(righttopnode,2)/2;

% point loads [GN]
BC.pointLoad = [];

% body force [GN/m3]
BC.b = @(x)[];  

%% Neumann BCs - fluid
% distributed flux [m/s]
% impervious at bottom, left, and right
BC.fluxNodes = [MeshP.left_dof; MeshP.right_dof; MeshP.bottom_dof; nodesP];
BC.fluxValue = zeros(length(BC.fluxNodes),1);

% point flux [m3/s]
BC.pointFlux = [];

% flux source [m3/s/m3]
BC.s = @(x)[]; 

%% Porosity BCs
if ~Control.Biotmodel
    BC.fixed_n = [];
    BC.free_n = setdiff(MeshN.DOF, BC.fixed_n);
    BC.fixed_n_value = zeros(length(BC.fixed_n),1);
end

%% Quadrature order
Control.nqU = 2;
Control.nqP = 2;

%% Problem type
% 1 = quasi-steady/transient problem (no acceleration and pressure change)
% 0 = dynamic problem (acceleration/intertia terms included)
Control.steady = 1;

% tag used for computing analytical solution
% 1 = uncoupled problem (elasticity, heat transfer, etc)
% 0 = coupled problem (Biot, Spanos model)
Control.uncoupled = 0; 

%% Solution parameters
Control.dt = 1e-2;  % time step
Control.tend = 10;   % final simulation time

Control.beta = 1; % beta-method time discretization -- beta = 1 Backward Euler; beta = 0.5 Crank-Nicolson

Control.plotu = 5*2; % dof y of node 5 (x = 5m, y = 10m)
Control.plotp = 5*2; % dof y of node 5 (x = 5m, y = 10m)

% plot analytical solution (valid for 1D problems with Material.Minv == 0)
Control.plotansol = 0; % 1 = true; 0 = false

% solve in the frequency domain
Control.freqDomain = 0;  % 1 = true; 0 = false

end