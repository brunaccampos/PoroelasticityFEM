function [Material, MeshU, MeshP, MeshN, BC, Control] = ManufacturedSolution1Dupv_Biot(~, progress_on,~,~)
% ------------------------------------------------------------------------
% Manufactured solution in 1D
% us = sin(x*t)/1000
% vf = cos(x*t)/1000
% p = sin(x*t)*cos(x*t)/1000
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
Control.PMmodel = 'Dyn6_Biot_UPV';

%% Material properties - Boone (1990) with modifications
% elasticity modulus [GPa]
Material.E = 14.4;
% Poisson's ratio
Material.nu = 0.2;
% intrinsic permeability [m2]
Material.k = 1.88e-13;
% dynamic viscosity [GPa s]
Material.mu = 1e-12;
% porous media permeability [m2/GPa s]
Material.kf = Material.k/Material.mu;
% Biot's coefficient
Material.alpha = 1;
% fluid bulk modulus [GPa]
Material.Kf = 3.3;
% solid bulk modulus [GPa]
Material.Ks = 36;
% material porosity
Material.eta0 = 0.19;
% 1/Q (related to storage coefficient)
Material.Minv = 0;
% fluid bulk viscosity [GPa s]
Material.xif = 2.8e-12; % (Quiroga-Goode, 2005)
% fluid density [10^9 kg/m3]
Material.rho_f = 1000e-9;
% solid density [10^9 kg/m3]
Material.rho_s = 2600e-9;
% average density of the medium
Material.rho = Material.eta0*Material.rho_f + (1-Material.eta0)*Material.rho_s;

% additional coefficients for analytical result
% Lame constant [GPa]
Material.lambda = Material.E * Material.nu/((1+Material.nu)*(1-2*Material.nu));
% gravitational acceleration [m/s2]
Material.g = 9.81;
% hydraulic conductivity [m/s]
Material.kh = Material.kf * Material.rho_f * Material.g;

% thickness
% 1D: cross sectional area [m2]
% 2D: out of plane thickness [m]
Material.t = 1;

% constititive law - 'PlaneStress' or 'PlaneStrain'
% Note: use 'PlaneStrain' for 1D or 2D poroelasticity
Material.constLaw = 'PlaneStress';

%% Spanos material parameters
% porosity effective pressure coefficient (Spanos, 1989)
% n = 0; % lower limit
n = 1; % return to Biot
% n = Material.Ks/Material.Kf; % upper limit

% modified storage coefficient (Muller, 2019)
Mstarinv = Material.Minv - (1-n)*(Material.alpha - Material.eta0)/Material.Ks; 
Mstar = 1/Mstarinv;

Material.deltaF = (Material.alpha - Material.eta0) * Material.eta0 * Mstar * n / Material.Ks;
Material.deltaS = (Material.alpha - Material.eta0) * Material.eta0 * Mstar /Material.Kf;

%% Mesh parameters
if progress_on
    disp([num2str(toc),': Building Mesh...']);
end

% location of initial node [m] [x0;y0;z0]
coord0 = [0;0;0];
% number of space dimensions
nsd = 1;
% size of domain [m] [Lx;Ly;Lz]
L = 4;
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

%% Initial conditions
% udot
BC.initU = ones(MeshU.nDOF,1)./1000;
% ufdot
BC.initUfdot = ones(MeshU.nDOF,1)./1000;

%% Dirichlet BCs - solid
% displacement prescribed on the left and right
BC.fixed_u = [MeshU.left_nodes; MeshU.right_nodes];
BC.fixed_u_value = @(t) [0; sin(L*t)]./1000;
% free displacement nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid velocity
% displacement prescribed on the left and right
BC.fixed_ufdot = [MeshU.left_nodes; MeshU.right_nodes];
BC.fixed_ufdot_value = @(t) [1; cos(L*t)]./1000;
% free displacement nodes
BC.free_ufdot = setdiff(MeshU.DOF, BC.fixed_ufdot);

%% Dirichlet BCs - fluid pressure
% pressure prescribed on the left and right
BC.fixed_p = [MeshP.left_nodes; MeshP.right_nodes];
BC.fixed_p_value = @(t) [0; sin(L*t)*cos(L*t)]./1000;
% free pressure nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% point load [GN]
BC.pointLoad = [];

% distributed load [GN/m2]
BC.tractionNodes = [];

% body force [GN/m3]
BC.bs = @(x,t) (Material.E * t^2 * sin(x*t) + (Material.alpha-Material.eta0)*t*(cos(x*t)^2-sin(x*t)^2) -...
    (1-Material.eta0)*Material.rho_s*x^2*sin(x*t) - Material.mu*Material.eta0^2/Material.k*(1-x)*cos(x*t))./1000;

BC.bf = @(x,t) (Material.eta0*(cos(x*t)^2-sin(x*t)^2) - Material.eta0*Material.rho_f*x*sin(x*t) + ...
    Material.mu*Material.eta0^2/Material.k*(1-x)*cos(x*t))./1000;

%% Neumann BCs - fluid
% point flux [m/s]
BC.pointFlux = [];

% distributed flux [m3/s]
BC.fluxNodes = [];

% flux source [m3/s/m3]
BC.s = @(x,t) (Material.Minv*x*(cos(x*t)^2-sin(x*t)^2) + (Material.alpha-Material.eta0)*(cos(x*t)-x*t*sin(x*t)) - ...
    Material.eta0*t*sin(x*t))./1000;

%% Quadrature order
Control.nqU = 3;
Control.nqP = 3;

%% Frequency domain
Control.freqDomain = 0;  % 1 = true; 0 = false

%% Analytical solution
% 1 = uncoupled problem (elasticity, heat transfer, etc)
% 0 = coupled problem (Biot, Spanos model)
Control.uncoupled = 1;

% plot analytical solution (valid for 1D problems with Material.Minv == 0)
Control.plotansol = 1; % 1 = true; 0 = false

% type of analytical solution to compute
% 'getAnSol_uncoupled' = uncoupled problem (elasticity, heat transfer, etc)
% 'getAnSol_coupledComp' = coupled porous media problem, compressible
% materials
% 'getAnSol_coupledIncomp' = coupled porous media problem, incompressible
% materials (1/M=0)
Control.ansol_type = 'getAnSol_uncoupled_upv';

% solution in u
Control.uan_symb = @(x,t) sin(x*t)./1000;
Control.u_an = @(t) Control.uan_symb(MeshU.coords,t);

Control.ufdotan_symb = @(x,t) cos(x*t)./1000;
Control.ufdot_an = @(t) Control.ufdotan_symb(MeshU.coords,t);

% solution in p
Control.pan_symb = @(x,t) (sin(x*t).*cos(x*t))./1000;
Control.p_an = @(t) Control.pan_symb(MeshP.coords,t);

%% Time step controls
Control.dt = 1e-3;  % time step [s]
Control.tend = 1;   % final simulation time [s]

% Newmark method
Control.beta = 0.7;
Control.gamma = 0.7;
Control.theta = 0.7;
Control.lambda = 0.7;

%% Plot data
% DOF to plot graphs
Control.plotu = round(length(MeshU.coords)/2);
Control.plotp = round(length(MeshP.coords)/2);

% Plot in a row
Control.fixedDepthPlotON = 0; % 0: false, 1: true

% Nodes to plot in a row (all nodes for 1D case)
Control.ploturow = MeshU.DOF;
Control.plotprow = MeshP.DOF;

end