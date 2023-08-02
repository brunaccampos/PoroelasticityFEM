function [Material, MeshU, MeshP, MeshN, BC, Control] = Column2D_Steady_Zheng(config_dir, progress_on)
% Column Consolidation 1D simulation
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
% column top at x=L, column bottom at x=0
% ------------------------------------------------------------------------

%% Poroelasticity model
% Options:  Tr1_Biot_UP -------- Biot model (u-p), transient
%           Tr2_Spanos_UPN ----- Spanos model (u-p-n), transient
%           Tr3_Spanos_UP ------ Spanos model (u-p), dynamic, implicit
%                                   porosity perturbation equation
%           Dyn1_Biot_UP -------- Biot model (u-p), dynamic
%           Dyn2_Spanos_UPN ----- Spanos model (u-p-n), dynamic
%           Dyn3_Spanos_UP ------ Spanos model (u-p), dynamic, implicit
%                                   porosity perturbation equation
%           Dyn4_Biot_UPU ------- Biot model (u-p-U), dynamic
%           Dyn5_Spanos_UPU ----- Spanos model (u-p-U), dynamic, implicit
%                                   porosity perturbation equation
Control.PMmodel = 'Tr1_Biot_UP';

%% Material properties - Berea Sandstone (Zheng, 2003)
% elasticity modulus [GPa]
Material.E = 14.4;
% Poisson's ratio
Material.nu = 0.2;
% intrinsic permeability [m2]
Material.k = 2e-13;
% dynamic viscosity [GPa s]
Material.mu = 1e-12;
% porous media permeability [m2/GPa s]
Material.kf = Material.k/Material.mu;
% Biot's coefficient
Material.alpha = 0.79;
% fluid bulk modulus [GPa]
Material.Kf = 3.3;
% solid bulk modulus [GPa]
Material.Ks = 36;
% material porosity
Material.n = 0.2;
% 1/Q (related to storage coefficient)
Material.Minv = (Material.alpha - Material.n)/Material.Ks + Material.n/Material.Kf;
% fluid bulk viscosity [GPa s]
Material.xif = 2.8e-12; % (Quiroga-Goode, 2005)

% lumped mass matrix - 0: false, 1: true
Material.lumpedMass = 0;

% thickness 
% 1D: cross sectional area [m2]
% 2D: out of plane thickness [m]
Material.t = 1;

% constititive law - 'PlaneStress' or 'PlaneStrain'
% Note: use 'PlaneStrain' for 1D or 2D poroelasticity
Material.constLaw = 'PlaneStrain';

%% Spanos material parameters
% porosity effective pressure coefficient (Spanos, 1989)
% n = 0; % lower limit
n = 1; % return to Biot
% n = Material.Ks/Material.Kf; % upper limit

% modified storage coefficient (Muller, 2019)
Mstarinv = Material.Minv - (1-n)*(Material.alpha - Material.n)/Material.Ks; 
Mstar = 1/Mstarinv;

Material.deltaF = (Material.alpha - Material.n) * Material.n * Mstar * n / Material.Ks;
Material.deltaS = (Material.alpha - Material.n) * Material.n * Mstar /Material.Kf;

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
        L = 3;
        %%%% solid displacement field
        typeU = 'L3';
        MeshU = Build1DMesh(nsd, ne, L, typeU);
        %%%% fluid pressure field
        typeP = 'L2';
        MeshP = Build1DMesh(nsd, ne, L, typeP);
        %%%% porosity field
        if contains(Control.PMmodel, 'UPN')
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
        meshFileNameU = 'Mesh Files\Plate_Zheng_Q9.msh';
        MeshU = BuildMesh_GMSH(meshFileNameU, fieldU, nsd, config_dir, progress_on);
        %%%% pressure field
        fieldP = 'p';
        meshFileNameP = 'Mesh Files\Plate_Zheng_Q4.msh';
        MeshP = BuildMesh_GMSH(meshFileNameP, fieldP, nsd, config_dir, progress_on);
        %%%% porosity field
        if contains(Control.PMmodel, 'UPN')
            fieldN = 'n';
            meshFileNameN = 'Mesh Files\Plate_Zheng_Q4.msh';
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
BC.traction = -4e-3;
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

%% Porosity BCs
if contains(Control.PMmodel, 'UPN')
    BC.fixed_n = [];
    BC.free_n = setdiff(MeshN.DOF, BC.fixed_n);
    BC.fixed_n_value = zeros(length(BC.fixed_n),1);
end

%% Quadrature order
Control.nqU = 2;
Control.nqP = 2;

%% Frequency domain
Control.freqDomain = 0;  % 1 = true; 0 = false

%% Analytical solution
% 1 = uncoupled problem (elasticity, heat transfer, etc)
% 0 = coupled problem (Biot, Spanos model)
Control.uncoupled = 0; 

% plot analytical solution (valid for 1D problems with Material.Minv == 0)
Control.plotansol = 0; % 1 = true; 0 = false

%% Time step controls
Control.dt = 1e-2;  % time step
Control.tend = 10;   % final simulation time

Control.beta = 1; % beta-method time discretization -- beta = 1 Backward Euler; beta = 0.5 Crank-Nicolson

%% Plot graphs
% DOF to plot graphs
Control.plotu = 388; % dof y of node 194 (x = 1.11m, y = 1.58m)
Control.plotp = 138; % dof of node 138 (x = 1.11m, y = 1.58m)

end