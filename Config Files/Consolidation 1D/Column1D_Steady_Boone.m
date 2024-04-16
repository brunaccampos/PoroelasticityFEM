function [Material, MeshU, MeshP, MeshN, BC, Control] = Column1D_Steady_Boone(~, progress_on,~,~)
% Column Consolidation 1D simulation
% Configuration File
% ------------------------------------------------------------------------
% Based on Korsawe (2006) model for transient/quasi-steady case
% ------------------------------------------------------------------------
% Assumptions/conventions:
% - stress is positive for tension
% - boundary condition for force is based on total stress
% - no acceleration terms for solid or fluid
% - solid velocity is neglected
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

%% Material properties - Boone (1990)
% shear modulus [GPa]
Material.mu = 6;
% Poisson's ratio
Material.nu = 0.2;
% elasticity modulus [GPa]
Material.E = 2 * Material.mu * (1 + Material.nu);
% porous media permeability [m2/GPa s]
Material.kf = 2e-2;
% dynamic viscosity [GPa s]
Material.muf = 1e-12;
% intrinsic permeability [m2]
Material.k = Material.kf * Material.muf;
% Biot's coefficient
Material.alpha = 1;
% fluid bulk modulus [GPa]
Material.Kf = 3;
% solid bulk modulus [GPa]
Material.Ks = 36;
% fluid bulk viscosity [GPa s]
Material.xif = 2.8e-12; % (Quiroga-Goode, 2005)
% material porosity
Material.eta0 = 0.19;
% 1/Q (related to storage coefficient)
Material.Minv = (Material.alpha - Material.eta0)/Material.Ks + Material.eta0/Material.Kf;

% additional coefficients for analytical result
% Lame constant [GPa]
Material.lambda = Material.E * Material.nu/((1+Material.nu)*(1-2*Material.nu));
% gravitational acceleration [m/s2]
Material.g = 9.81;
% fluid density [10^9 kg/m3]
Material.rhof = 1000e-9;
% hydraulic conductivity [m/s]
Material.kh = Material.kf * Material.rhof * Material.g;

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
Mstarinv = Material.Minv - (1-n)*(Material.alpha - Material.eta0)/Material.Ks; 
Mstar = 1/Mstarinv;

% porosity equation coefficients
Material.deltaf = (Material.alpha - Material.eta0) * Material.eta0 * Mstar * n / Material.Ks;
Material.deltas = (Material.alpha - Material.eta0) * Material.eta0 * Mstar / Material.Kf;

%% Mesh parameters
if progress_on
    disp([num2str(toc),': Building Mesh...']);
end

% location of initial node [m] [x0;y0;z0]
coord0 = [0;0;0];
% number of space dimensions
nsd = 1;
% size of domain [m] [Lx;Ly;Lz]
L = 6;
% number of elements in each direction [nex; ney; nez]
ne = 100;

%%%% displacement mesh
% element type ('Q4')
typeU = 'L3';
% variable field ('u', 'p', 'n')
fieldU = 'u';
MeshU = BuildMesh_structured(nsd, coord0, L, ne, typeU, fieldU, progress_on);

%%%% pressure mesh
% element type ('Q4')
typeP = 'L2';
% variable field ('u', 'p', 'n')
fieldP = 'p';
MeshP = BuildMesh_structured(nsd, coord0, L, ne, typeP, fieldP, progress_on);

%%%% porosity mesh
if contains(Control.PMmodel, 'UPN')
    % element type ('Q4')
    typeN = 'L2';
    % variable field ('u', 'p', 'n')
    fieldN = 'n';
    MeshN = BuildMesh_structured(nsd, coord0, L, ne, typeN, fieldN, progress_on);
else
    MeshN = [];
end

%% Dirichlet BCs - solid
% displacement u=0 at the bottom
BC.fixed_u = MeshU.right_nodes;
BC.fixed_u_value = @(t) zeros(length(BC.fixed_u),1);
% free displacement nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid
% pressure p=0 at the top
BC.fixed_p = MeshP.left_nodes;
BC.fixed_p_value = @(t) zeros(length(BC.fixed_p),1);
% free pressure nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% point load [GN]
BC.pointLoadValue = 1e-6;
BC.pointLoadNodes = MeshU.left_nodes;
BC.pointLoad = zeros(MeshU.nDOF,1);
BC.pointLoad(BC.pointLoadNodes) = BC.pointLoadValue;

% distributed load [GN/m2]
BC.tractionNodes = [];

% body force [GN/m3]
BC.b = @(x,t)[];  

%% Neumann BCs - fluid
% point flux [m/s]
BC.pointFluxValue = 0;
BC.pointFluxNodes = [];
BC.pointFlux = zeros(MeshP.nDOF,1);
BC.pointFlux(BC.pointFluxNodes) = BC.pointFluxValue;

% distributed flux [m3/s]
BC.fluxNodes = [];

% flux source [m3/s/m3]
BC.s = @(x,t)[]; 

%% Porosity BCs
if contains(Control.PMmodel, 'UPN')
    BC.fixed_n = [];
    BC.free_n = setdiff(MeshN.DOF, BC.fixed_n);
    BC.fixed_n_value = zeros(length(BC.fixed_n),1);
end

%% Quadrature order
Control.nqU = 2;
Control.nqP = 1;

%% Frequency domain
Control.freqDomain = 0;  % 1 = true; 0 = false

%% Solution parameters
% tag used for computing analytical solution
% 1 = uncoupled problem (elasticity, heat transfer, etc)
% 0 = coupled problem (Biot, Spanos model)
Control.uncoupled = 0; 

% plot analytical solution (valid for 1D problems)
Control.plotansol = 1; % 1 = true; 0 = false

% type of analytical solution to compute
% 'getAnSol_uncoupled' = uncoupled problem (elasticity, heat transfer, etc)
% 'getAnSol_coupledComp' = coupled porous media problem, compressible
% materials
% 'getAnSol_coupledIncomp' = coupled porous media problem, incompressible
% materials (1/M=0)
Control.ansol_type = 'getAnSol_coupledComp';

Control.u_an = @(t) zeros(MeshU.nDOF,1);
Control.p_an = @(t) zeros(MeshP.nDOF,1);

%% Time step controls
Control.dt = 1e-2;  % time step [s]
Control.tend = 10;   % final simulation time [s]

Control.beta = 1; % beta-method time discretization -- beta = 1 Backward Euler; beta = 0.5 Crank-Nicolson

%% Plot data
% DOF to plot graphs
Control.plotu = round(length(MeshU.coords)/2);
Control.plotp = round(length(MeshP.coords)/2);

% Plot in a row
Control.fixedDepthPlotON = 1; % 0: false, 1: true

% Nodes to plot in a row (all nodes for 1D case)
Control.ploturow = MeshU.DOF;
Control.plotprow = MeshP.DOF;

end