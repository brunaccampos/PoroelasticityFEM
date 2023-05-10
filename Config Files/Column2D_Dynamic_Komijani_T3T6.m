function [Material, MeshU, MeshP, MeshN, BC, Control] = Column2D_Dynamic_Komijani_T3T6(config_dir, progress_on)
% Column Consolidation 2D simulation
% Configuration File
% Based on Zienkiewicz (1982) model
% ------------------------------------------------------------------------
% Assumptions/conventions:
% - stress is positive for tension
% - boundary condition for force is based on total stress
% - only solid acceleration is considered (undrained condition; no motions
% of the fluid relative to the solid skeleton can occur)
% - solid grains and fluid are incompressible
% ------------------------------------------------------------------------
% column top at x=L, column bottom at x=0
% ------------------------------------------------------------------------

%% Poroelasticity model
% Options:  Transient_Biot ----- Biot model (u-p), transient
%           Transient_Spanos --- Spanos model (u-p-n), transient
%           Transient_BiotPoro - Biot model (u-p), dynamic, implicit
%                                   porosity perturbation equation
%           Dynamic_Biot ------- Biot model (u-p), dynamic
%           Dynamic_Spanos ----- Spanos model (u-p-n), dynamic
%           Dynamic_BiotPoro --- Biot model (u-p), dynamic, implicit
%                                   porosity perturbation equation
Control.Biotmodel = 'Dynamic_Biot';

%% Material properties - Komijani (2019)
% elasticity modulus [GPa]
Material.E = 14.516e-3;
% Poisson's ratio
Material.nu = 0.3;
% porous media permeability [m2/GPa s]
Material.kf = 1.0194e3;
% dynamic viscosity [GPa s]
Material.mu = 1e-12;
% intrinsic permeability [m2]
Material.k = Material.kf * Material.mu;
% fluid bulk modulus [GPa]
Material.Kf = 2.1;
% solid bulk modulus [GPa]
Material.Ks = 1e11;
% material porosity
Material.n = 0.3;
% Biot's coefficient
Material.alpha = 1;
% fluid density [10^9 kg/m3]
Material.rho_f = 1000e-9;
% solid density [10^9 kg/m3]
Material.rho_s = 2000e-9;
% average density of the medium
Material.rho = Material.n*Material.rho_f + (1-Material.n)*Material.rho_s;
% 1/Q (related to storage coefficient)
Material.Minv = (Material.alpha - Material.n)/Material.Ks + Material.n/Material.Kf;

% thickness 
% 1D: cross sectional area [m2]
% 2D: out of plane thickness [m]
Material.t = 1;

% constititive law - 'PlaneStress' or 'PlaneStrain'
% Note: use 'PlaneStrain' for 1D or 2D poroelasticity
Material.constLaw = 'PlaneStrain';

% lumped mass matrix - 0: false, 1: true
Material.lumpedMass = 0;

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
        if contains(Control.Biotmodel, 'Spanos')
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
        meshFileNameU = 'Mesh Files\Column2DT6.msh';
        MeshU = BuildMesh_GMSH(meshFileNameU, fieldU, nsd, config_dir, progress_on);
        %%%% pressure field
        fieldP = 'p';
        meshFileNameP = 'Mesh Files\Column2DT3.msh';
        MeshP = BuildMesh_GMSH(meshFileNameP, fieldP, nsd, config_dir, progress_on);
        %%%% porosity field
        if contains(Control.Biotmodel, 'Spanos')
            fieldN = 'n';
            meshFileNameN = 'Mesh Files\Column2DT3.msh';
            MeshN = BuildMesh_GMSH(meshFileNameN, fieldN, nsd, config_dir, progress_on);
        else
            MeshN = [];
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
BC.fixed_p = [MeshP.top_dof];
% fixed DOF values
BC.fixed_p_value = zeros(length(BC.fixed_p),1);
% free nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% traction interpolation (needed for traction applied in wells); 1 - true, 0 - false
BC.tractionInterp = 0;
% prescribed traction [GN/m2]
BC.traction = -3000e-9;
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

% point loads [GN]
BC.pointLoad = [];

% body force [GN/m3]
BC.b = @(x)[];  

%% Neumann BCs - fluid
% distributed flux [m/s]
% impervious at bottom, left, and right
BC.fluxNodes = [MeshP.left_dof; MeshP.right_dof; MeshP.bottom_dof];
BC.fluxValue = zeros(length(BC.fluxNodes),1);

% point flux [m3/s]
BC.pointFlux = [];

% flux source [m3/s/m3]
BC.s = @(x)[]; 

%% Quadrature order
Control.nqU = 2;
Control.nqP = 2;

%% Solution parameters
% tag used for computing analytical solution
% 1 = uncoupled problem (elasticity, heat transfer, etc)
% 0 = coupled problem (Biot, Spanos model)
Control.uncoupled = 0; 

Control.dt = 1e-2;  % time step
Control.tend = 10;   % final simulation time
Control.tol = 1e-3; % tolerance for NR method
Control.max_it = 100; % maximum of iterations

Control.plotu = 200; % dof y of node 100 (x = 0.05m, y = 5m)
Control.plotp = 56; % dof of node 56 (x = 0.05m, y = 5m)

% Control.plotu = 32; % dof y of node 16 (x = 0.1m, y = 5m)
% Control.plotp = 14; % dof of node 14 (x = 0.1m, y = 5m)

% plot analytical solution (valid for 1D problems with Material.Minv == 0)
Control.plotansol = 0; % 1 = true; 0 = false

% solve in the frequency domain
Control.freqDomain = 0;  % 1 = true; 0 = false

%% Time discretization parameters
% Newmark method
Control.beta = 0.7;
Control.gamma = 0.7;
Control.theta = 0.7;

end