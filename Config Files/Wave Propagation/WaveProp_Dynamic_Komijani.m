function [Material, MeshU, MeshP, MeshN, BC, Control] = WaveProp_Dynamic_Komijani(config_dir, progress_on,~,~)
% Wave propagation in 2D
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

%% Material properties - Quiroga-Goode (2005)
% elasticity modulus [Pa]
Material.E = 200e9;
% Poisson's ratio
Material.nu = 0.3;
% average density of the medium [kg/m3]
Material.rho = 8000;

% porous media permeability [m2/GPa s]
Material.kf = 0;
% Biot's coefficient
Material.alpha = 0;
% 1/Q (related to storage coefficient)
Material.Minv = 0;
% fluid density [10^9 kg/m3]
Material.rhof = 0;

% lumped mass matrix - 0: false, 1: true
Material.lumpedMass = 1;

% thickness 
% 1D: cross sectional area [m2]
% 2D: out of plane thickness [m]
Material.t = 1;

% constititive law - 'PlaneStress' or 'PlaneStrain'
% Note: use 'PlaneStrain' for 1D or 2D poroelasticity
Material.constLaw = 'PlaneStress';

%% Mesh parameters
if progress_on
    disp([num2str(toc),': Building Mesh...']);
end

% GMSH file Version 2 ASCII
% number of space dimensions
nsd = 2;
%%%% displacement field
fieldU = 'u';
meshFileNameU = 'Mesh Files\WavePropPlate_100x150cmQ4.msh';
MeshU = BuildMesh_GMSH(meshFileNameU, fieldU, nsd, config_dir, progress_on);
%%%% pressure field
fieldP = 'p';
meshFileNameP = 'Mesh Files\WavePropPlate_100x150cmQ4.msh';
MeshP = BuildMesh_GMSH(meshFileNameP, fieldP, nsd, config_dir, progress_on);
%%%% porosity field
if contains(Control.PMmodel, 'UPN')
    fieldN = 'n';
    meshFileNameN = 'Mesh Files\WavePropPlate_100x150cmQ4.msh';
    MeshN = BuildMesh_GMSH(meshFileNameN, fieldN, nsd, config_dir, progress_on);
else
    MeshN = [];
end

%% Dirichlet BCs - solid
% displacement fixed at central node in y
BC.fixed_u = MeshU.left_dofx;
t0 = 1e-5;
% fixed DOF values
BC.fixed_u_value = @(t) ones(length(BC.fixed_u),1) * 0.03*sin(pi*(t<=(2*t0))/2/t0);
% free displacement nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid
%   pressure p=0 at the top
BC.fixed_p = 1:MeshP.nDOF;
% fixed DOF values
BC.fixed_p_value = @(t) zeros(length(BC.fixed_p),1);
% free pressure nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% point load [GN]
BC.pointLoad = [];

% distributed load [GN/m2]
BC.tractionNodes = [];

% body force [GN/m3]
BC.b = @(x,t)[];  

%% Neumann BCs - fluid
% point flux [m/s]
BC.pointFlux = [];

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
Control.nqP = 2;

%% Frequency domain
Control.freqDomain = 0;  % 1 = true; 0 = false

%% Analytical solution
% plot analytical solution (valid for 1D problems with Material.Minv == 0)
Control.plotansol = 0; % 1 = true; 0 = false

%% Time step controls
Control.dt = 188e-8;  % time step
Control.tend = 188e-6;   % final simulation time

% Newmark method
Control.beta = 0.7;
Control.gamma = 0.7;
Control.theta = 0.7;
Control.lambda = 0.7;

%% Plot data
% DOF to plot graphs
Control.plotu = 611*2; % dof y of node 242 (x = 7.5m, y = 7.5m)
Control.plotp = 491; % dof y of node 177 (x = 7.5m, y = 7.5m)

end