function [Material, MeshU, MeshP, MeshN, BC, Control] = ManufacturedSolution1D_Biot(config_dir, progress_on,~,~)
% ------------------------------------------------------------------------
% Manufactured solution in 1D
% u = sin(xt)
% p = cos(xt)
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
Material.G = 1e3;
% Poisson's ratio
Material.nu = 0.3;
% elasticity modulus [GPa]
Material.E = 2 * Material.G * (1 + Material.nu);
% intrinsic permeability [m2]
Material.k = 0.1;
% dynamic viscosity [GPa s]
Material.mu = 1e-3;
% porous media permeability [m2/GPa s]
Material.kf = Material.k/Material.mu;
% Biot's coefficient
Material.alpha = 1;
%  [GPa]
Material.Ku = 1.5;
%  [GPa]
Material.B = 0.9;
% 1/Q (related to storage coefficient)
% Material.Minv = Material.alpha/(Material.Ku*Material.B);
Material.Minv = 0.75;
% additional coefficients for analytical result
% Lame constant [GPa]
Material.lambda = Material.E * Material.nu/((1+Material.nu)*(1-2*Material.nu));
% gravitational acceleration [m/s2]
Material.g = 9.81;
% fluid density [10^9 kg/m3]
Material.rho_f = 1000e-9;
% hydraulic conductivity [m/s]
Material.kh = Material.kf * Material.rho_f * Material.g;

% thickness
% 1D: cross sectional area [m2]
% 2D: out of plane thickness [m]
Material.t = 1;

% constititive law - 'PlaneStress' or 'PlaneStrain'
% Note: use 'PlaneStrain' for 1D or 2D poroelasticity
Material.constLaw = 'PlaneStrain';

%% Mesh parameters
if progress_on
    disp([num2str(toc),': Building Mesh...']);
end

% location of initial node [m] [x0;y0;z0]
coord0 = [0;0;0];
% number of space dimensions
nsd = 1;
% size of domain [m] [Lx;Ly;Lz]
L = 5;
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
% pressure
BC.initP = zeros(MeshP.nDOF,1);
% pressure p=1 at the top
BC.initP(MeshP.left_nodes,1) = 1;

%% Dirichlet BCs - solid
% displacement u=0 at the top
BC.fixed_u = (MeshU.left_nodes);
BC.fixed_u_value = @(t) zeros(length(BC.fixed_u),1);
% free displacement nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid
% pressure p=1 at the top
BC.fixed_p = (MeshP.left_nodes);
BC.fixed_p_value = @(t) ones(length(BC.fixed_p),1);
% free pressure nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% point load [GN]
L = max(MeshU.coords);
BC.pointLoadNodes = MeshU.right_nodes;
BC.pointLoad = @(t) [zeros(MeshU.nDOF-1,1); (2*Material.G + Material.lambda)* t * cos(L*t) - Material.alpha * cos(L*t)];

% distributed load [GN/m2]
BC.tractionNodes = [];

% body force [GN/m3]
BC.b = @(x,t) (2*Material.G + Material.lambda) * t^2 * sin (x*t) - Material.alpha * t * sin(x*t);

%% Neumann BCs - fluid
% point flux [m/s]
BC.pointFluxNodes = MeshP.right_nodes;
BC.pointFlux = @(t) [zeros(MeshP.nDOF-1,1); t * Material.kf * sin(L*t)];

% distributed flux [m3/s]
BC.fluxNodes = [];

% flux source [m3/s/m3]
BC.s = @(x,t) Material.alpha * (cos(x*t) - t * x * sin(x*t)) - Material.Minv * x * sin(x*t) + Material.kf * t^2 * cos(x*t);

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
Control.ansol_type = 'getAnSol_uncoupled';

% solution in u
Control.uan_symb = @(x,t) sin(x*t);
Control.u_an = @(t) Control.uan_symb(MeshU.coords,t);

% solution in p
Control.pan_symb = @(x,t) cos(x*t);
Control.p_an = @(t) Control.pan_symb(MeshP.coords,t);

%% Time step controls
Control.dt = 1e-2;  % time step [s]
Control.tend = 1;   % final simulation time [s]

Control.beta = 1; % beta-method time discretization -- beta = 1 Backward Euler; beta = 0.5 Crank-Nicolson

%% Plot data
% DOF to plot graphs
Control.plotu = round(length(MeshU.coords)/2);
Control.plotp = round(length(MeshP.coords)/2);

end