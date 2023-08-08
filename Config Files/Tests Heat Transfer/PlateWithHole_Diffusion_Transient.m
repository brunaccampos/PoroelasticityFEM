function [Material, MeshU, MeshP, MeshN, BC, Control] = PlateWithHole_Diffusion_Transient(config_dir, progress_on)
% 2D diffusion problem - transient case
% Configuration file
% ------------------------------------------------------------------------
% Adapted from: https://github.com/GCMLab
% ------------------------------------------------------------------------

%% Material properties
% diffusion coefficient [m2/s] 
Material.kf = 1e-3;
% 1/Q (related to storage coefficient)
Material.Minv = 1;

% material density [kg/m3]
Material.rho = 0;
% fluid density [kg/m3]
Material.rho_f = 0;
% elasticity modulus [Pa]
Material.E = 0;
% Poisson's ratio
Material.nu = 0;
% Biot's coefficient
Material.alpha = 0;
% poroelasticity model
Control.PMmodel = 'Tr1_Biot_UP';

% thickness 
% 1D: cross sectional area [m2]
% 2D: out of plane thickness [m]
Material.t = 1;

% constititive law - 'PlaneStress' or 'PlaneStrain'
% Note: use 'PlaneStrain' for 1D or 2D poroelasticity
Material.constLaw = 'PlaneStress';

%% Mesh Properties
if progress_on
    disp([num2str(toc),': Building Mesh...']);
end

% GMSH file Version 2 ASCII
% number of space dimensions
nsd = 2;
%%%% displacement field
fieldU = 'u';
meshFileNameU = 'Mesh Files\PlateWithHoleQ4.msh';
MeshU = BuildMesh_GMSH(meshFileNameU, fieldU, nsd, config_dir, progress_on);
%%%% pressure field
fieldP = 'p';
meshFileNameP = 'Mesh Files\PlateWithHoleQ4.msh';
MeshP = BuildMesh_GMSH(meshFileNameP, fieldP, nsd, config_dir, progress_on);
%%%% porosity field
if contains(Control.PMmodel, 'UPN')
    fieldN = 'n';
    meshFileNameN = 'Mesh Files\PlateWithHoleQ4.msh';
    MeshN = BuildMesh_GMSH(meshFileNameN, fieldN, nsd, config_dir, progress_on);
else
    MeshN = [];
end

%% Dirichlet BCs - solid
% column vector of prescribed displacement dof
BC.fixed_u = 1:MeshU.nDOF;
% prescribed displacement for each dof [u1; u2; ...] [m]
BC.fixed_u_value = @(t) zeros(length(BC.fixed_u),1);
% free nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid
T = 100;
BC.fixed_p = [1; 5; 6; 7; 8; 9; 10; 11; 12]; % nodes at the inner circle
% fixed DOF values
BC.fixed_p_value = @(t) T*ones(length(BC.fixed_p),1);
% free nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% distributed traction [N/m2]
BC.tractionNodes = 1:MeshU.nDOF;
BC.tractionForce = zeros(length(BC.tractionNodes),2);

% point loads [N]
BC.pointLoad = [];

% body force [N/m3]
BC.b = @(x,t)[];  

%% Neumann BCs - fluid
% distributed flux [m/s]
BC.fluxNodes = [];

% point flux [m3/s]
BC.pointFlux = [];

% flux source [N/m3]
BC.s = @(x,t)[]; 

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
Control.dt = 0.1;  % time step
Control.tend = 10;   % final simulation time

Control.beta = 1; % beta-method time discretization -- beta = 1 Backward Euler; beta = 0.5 Crank-Nicolson

% Newmark method
Control.beta = 0.7;
Control.gamma = 0.7;
Control.theta = 0.7;

%% Plot data
% DOF to plot graphs
Control.plotu = 2;
Control.plotp = 2;

end