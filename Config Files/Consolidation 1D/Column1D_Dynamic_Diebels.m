function [Material, MeshU, MeshP, MeshN, BC, Control] = Column1D_Dynamic_Diebels(~, progress_on,~,~)
% Column Consolidation 1D simulation
% Configuration File
% ------------------------------------------------------------------------
% Based on Zienkiewicz (1982) model for dynamic case
% ------------------------------------------------------------------------
% Assumptions/conventions:
% - stress is positive for tension
% - boundary condition for force is based on total stress
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
Control.PMmodel = 'Dyn1_Biot_UP';

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
Material.muf = 1e-12;
% intrinsic permeability [m2]
Material.k = Material.kf * Material.muf;
% material porosity
Material.eta0 = 0.33;
% Biot's coefficient
Material.alpha = 1;
% fluid density [10^9 kg/m3]
Material.rhof = 1000e-9;
% solid density [10^9 kg/m3]
Material.rhos = 2000e-9;
% average density of the medium
Material.rho = Material.eta0*Material.rhof + (1-Material.eta0)*Material.rhos;
% 1/Q (related to storage coefficient)
Material.Minv = 0;

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
L = 10;
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
BC.pointLoadValue = 3000e-9;
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
BC.pointFluxNodes = MeshP.right_nodes;
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
Control.nqU = 3;
Control.nqP = 3;

%% Frequency domain
Control.freqDomain = 0;  % 1 = true; 0 = false

%% Analytical solution
% plot analytical solution (valid for 1D problems)
Control.plotansol = 0; % 1 = true; 0 = false

%% Time step controls
Control.dt = 1e-3;  % time step
Control.tend = 0.4;   % final simulation time

% Newmark method
Control.beta = 0.7;
Control.gamma = 0.7;
Control.theta = 0.7;
Control.lambda = 0.7;

%% Plot data
% DOF to plot graphs
Control.plotu = find(MeshU.coords == 5); % x = 5m
Control.plotp = find(MeshP.coords == 5); % x = 5m

end