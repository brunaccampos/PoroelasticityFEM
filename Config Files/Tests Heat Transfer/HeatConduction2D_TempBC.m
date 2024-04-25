function [Material, MeshU, MeshP, MeshN, BC, Control] = HeatConduction2D_TempBC(config_dir, progress_on,~,~)
% Plate with hole 1/8 model: Heat transfer problem adapted from file Q4one8thModel
% Configuration file
% ------------------------------------------------------------------------

%% Material properties
% thermal conductance coefficient [W/m3]
Material.kf = 1;
% Poroelasticity model
Control.PMmodel = 'Tr1_Biot_UP';

% elasticity modulus [Pa]
Material.E = 0;
% Poisson's ratio
Material.nu = 0;
% Biot's coefficient
Material.alpha = 0;
% 1/Q (related to storage coefficient)
Material.Minv = 0;
% poroelasticity model
Control.Biotmodel = 1;

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

% GMSh file Version 2 ASCII
% number of space dimensions
nsd = 2;
%%%% displacement field
fieldU = 'u';
meshFileNameU = 'Mesh Files\UnitPlateQ4.msh';
MeshU = BuildMesh_GMSH(meshFileNameU, fieldU, nsd, config_dir, progress_on);
%%%% pressure field
fieldP = 'p';
meshFileNameP = 'Mesh Files\UnitPlateQ4.msh';
MeshP = BuildMesh_GMSH(meshFileNameP, fieldP, nsd, config_dir, progress_on);
%%%% porosity field
if contains(Control.PMmodel, 'UPN')
    fieldN = 'n';
    meshFileNameN = 'Mesh Files\UnitPlateQ4.msh';
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
BC.fixed_p1 = MeshP.top_nodes;
BC.fixed_p2 = [MeshP.right_nodes; MeshP.left_nodes; MeshP.bottom_nodes];
BC.fixed_p = [BC.fixed_p1; BC.fixed_p2];
% fixed DOF values
BC.fixed_p_value = @(t) [250*ones(length(BC.fixed_p1),1); 40*ones(length(BC.fixed_p2),1)];

% free nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% distributed traction [N/m2]
BC.tractionNodes = [];

% point loads [N]
BC.pointLoad = @(t)[];

% body force [N/m3]
BC.b = @(x,t)[];  

%% Neumann BCs - fluid
% distributed flux [m/s]
BC.fluxNodes = [];

% point flux [m3/s]
BC.pointFlux = @(t)[];

% flux source [m3/s/m3]
BC.s = @(x,t)[]; 

%% Quadrature order
Control.nqU = 2;
Control.nqP = 2;

%% Frequency domain
Control.freqDomain = 0;  % 1 = true; 0 = false

%% Analytical solution
% plot analytical solution (valid for 1D problems with Material.Minv == 0)
Control.plotansol = 0; % 1 = true; 0 = false

%% Time step controls
Control.dt = 1;  % time step
Control.tend = 1;   % final simulation time

Control.beta = 1; % beta-method time discretization -- beta = 1 Backward Euler; beta = 0.5 Crank-Nicolson

%% Plot data
% DOF to plot graphs
Control.plotu = 1;
Control.plotp = 1;

end