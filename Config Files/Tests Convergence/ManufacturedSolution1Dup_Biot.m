function [Material, MeshU, MeshP, MeshN, BC, Control] = ManufacturedSolution1Dup_Biot(~, progress_on,~,~)
% ------------------------------------------------------------------------
% Manufactured solution in 1D
% us = sin(F1*x + F2*t)
% uf = cos(F1*x + F2*t)
% p = sin(F1*x) + sin(F2*t)
% ------------------------------------------------------------------------
% Based on Korsawe (2006) model for transient/quasi-steady case
% ------------------------------------------------------------------------
% Assumptions/conventions:
% - stress is positive for tension
% - boundary condition for force is based on total stress
% - no acceleration terms for solid or fluid
% - solid velocity is neglected
% ------------------------------------------------------------------------
% Porous media theories
% - BT: Biot
% - dCS: de la Cruz and Spanos
% ------------------------------------------------------------------------
% Loading options
% - Tr: transient/quasi-steady
% - Dyn: dynamic (acceleration included)
% ------------------------------------------------------------------------
% Main variables
% u = solid displacement
% p = fluid pressure
% n = porosity
% U = fluid displacement
% v = fluid velocity
% w = relative fluid velocity
% ------------------------------------------------------------------------
% Model options
%
% Tr_BT_UP          Tr_dCS_UP           Tr_dCS_UPN 
%
% Dyn_BT_UP         Dyn_BT_UPU          Dyn_BT_UPV          Dyn_BT_UPW
%
% Dyn_dCS_UP        Dyn_dCS_UPU         Dyn_dCS_UPN         Dyn_dCS_UPW
% ------------------------------------------------------------------------

%% Poroelasticity model
Control.PMmodel = 'Dyn_BT_UP';

%% Material properties - Boone (1990)
% elasticity modulus [GPa]
Material.M(1).E = 14.4;
% Poisson's ratio
Material.M(1).nu = 0.2;
% intrinsic permeability [m2]
Material.M(1).k = 1.88e-13;
% dynamic viscosity [GPa s]
Material.M(1).muf = 1e-12;
% porous media permeability [m2/GPa s]
Material.M(1).kf = Material.M(1).k/Material.M(1).muf;
% Biot's coefficient
Material.M(1).alpha = 0.79;
% fluid bulk modulus [GPa]
Material.M(1).Kf = 3.3;
% solid bulk modulus [GPa]
Material.M(1).Ks = 36;
% material porosity
Material.M(1).eta0 = 0.19;
% 1/Q (related to storage coefficient)
Material.M(1).Minv = (Material.M(1).alpha - Material.M(1).eta0)/Material.M(1).Ks + Material.M(1).eta0/Material.M(1).Kf;
% fluid bulk viscosity [GPa s]
Material.M(1).xif = 2.8e-12; % (Quiroga-Goode, 2005)
% fluid density [10^9 kg/m3]
Material.M(1).rhof = 1000e-9;
% solid density [10^9 kg/m3]
Material.M(1).rhos = 2600e-9;
% average density of the medium
Material.M(1).rho = Material.M(1).eta0*Material.M(1).rhof + (1-Material.M(1).eta0)*Material.M(1).rhos;

% additional coefficients for analytical result
% Lame constant [GPa]
Material.M(1).lambda = Material.M(1).E * Material.M(1).nu/((1+Material.M(1).nu)*(1-2*Material.M(1).nu));
% gravitational acceleration [m/s2]
Material.M(1).g = 9.81;
% hydraulic conductivity [m/s]
Material.M(1).kh = Material.M(1).kf * Material.M(1).rhof * Material.M(1).g;

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
% n = Material.M(1).Ks/Material.M(1).Kf; % upper limit

% modified storage coefficient (Muller, 2019)
Mstarinv = Material.M(1).Minv - (1-n)*(Material.M(1).alpha - Material.M(1).eta0)/Material.M(1).Ks; 
Mstar = 1/Mstarinv;

Material.M(1).deltaf = (Material.M(1).alpha - Material.M(1).eta0) * Material.M(1).eta0 * Mstar * n / Material.M(1).Ks;
Material.M(1).deltas = (Material.M(1).alpha - Material.M(1).eta0) * Material.M(1).eta0 * Mstar /Material.M(1).Kf;

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
ne = 192;

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

F1 = 1; % frequency 1 [Hz]
F2 = 2; % frequency 2 [Hz]

%% Initial conditions
% us
BC.initU = sin(F1.*MeshU.coords)./1000;
% p
BC.initP = sin(F1.*MeshP.coords)./1000;
% usdot
BC.initUdot = F2*cos(F1.*MeshU.coords)./1000;

%% Dirichlet BCs - solid
% displacement prescribed on the left and right
BC.fixed_u = [MeshU.left_nodes; MeshU.right_nodes];
BC.fixed_u_value = @(t) [sin(F2*t)*ones(length(MeshU.left_nodes),1); sin(F1*L+F2*t)*ones(length(MeshU.right_nodes),1)]./1000;
% free displacement nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid pressure
% pressure prescribed on the left and right
BC.fixed_p = [MeshP.left_nodes; MeshP.right_nodes];
BC.fixed_p_value = @(t) [sin(F2*t)*ones(length(MeshP.left_nodes),1); (sin(F1*L)+sin(F2*t))*ones(length(MeshP.right_nodes),1)]./1000;
% free pressure nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% point load [GN]
BC.pointLoad = @(t)[];

% distributed load [GN/m2]
BC.tractionNodes = [];

% body force [GN/m3]
BC.b = @(x,t) (Material.M(1).E*F1^2*sin(F1*x+F2*t) + Material.M(1).alpha*F1*cos(F1*x) - ...
    Material.M(1).rho*F2^2*sin(F1*x+F2*t))./1000;

%% Neumann BCs - fluid
% point flux [m/s]
BC.pointFlux = @(t)[];

% distributed flux [m3/s]
BC.fluxNodes = [];

% flux source [m3/s/m3]
BC.s = @(x,t) (Material.M(1).k/Material.M(1).muf*F1^2*sin(F1*x) + ...
    Material.M(1).k/Material.M(1).muf*Material.M(1).rhof*F2^2*F1*cos(F1*x+F2*t) - ...
    Material.M(1).alpha*F1*F2*sin(F1*x+F2*t) + Material.M(1).Minv*F2*cos(F2*t))./1000;


%% Porosity BCs
if contains(Control.PMmodel, 'UPN')
    BC.fixed_n = [];
    BC.free_n = setdiff(MeshN.DOF, BC.fixed_n);
    BC.fixed_n_value = zeros(length(BC.fixed_n),1);
end

%% Quadrature order
Control.nqU = 3;
Control.nqP = 3;

%% Frequency domain
Control.freqDomain = 0;  % 1 = true; 0 = false

%% Analytical solution
% plot analytical solution (valid for 1D problems)
Control.plotansol = 1; % 1 = true; 0 = false

% type of analytical solution to compute
% 'getAnSol_uncoupled' = uncoupled problem (elasticity, heat transfer, etc)
% 'getAnSol_coupledComp' = coupled porous media problem, compressible
% materials
% 'getAnSol_coupledIncomp' = coupled porous media problem, incompressible
% materials (1/M=0)
Control.ansol_type = 'getAnSol_uncoupled_UP';

% solution in u
Control.uan_symb = @(x,t) sin(F1*x + F2*t)./1000;
Control.u_an = @(t) Control.uan_symb(MeshU.coords,t);

% solution in p
Control.pan_symb = @(x,t) (sin(F1*x) + sin(F2*t))./1000;
Control.p_an = @(t) Control.pan_symb(MeshP.coords,t);

%% Time step controls
Control.dt = 1e-5;  % time step [s]
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