% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function [Material, MeshU, MeshP, MeshN, BC, Control] = PlateDiffusionDynamic(config_dir, progress_on,~,~)
% 2D diffusion problem - dynamic case

%% Material properties
% diffusion coefficient [ppm/mm2] 
Material.M(1).kf = 5e-3;
% 1/Q (related to storage coefficient)
Material.M(1).Minv = 1;
% poroelasticity model
Control.PMmodel = 'Dyn_BT_UP';

% material density [kg/m3]
Material.M(1).rho = 0;
% fluid density [kg/m3]
Material.M(1).rhof = 0;
% elasticity modulus [Pa]
Material.M(1).E = 0;
% Poisson's ratio
Material.M(1).nu = 0;
% Biot's coefficient
Material.M(1).alpha = 0;

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
meshFileNameU = 'Mesh Files\PlateDiffusion.msh';
MeshU = BuildMesh_GMSH(meshFileNameU, fieldU, nsd, config_dir, progress_on);
% type of material per element
MeshU.MatList = zeros(MeshU.ne, 1, 'int8');
% assign material type to elements
MeshU.MatList(:) = 1;

%%%% pressure field
fieldP = 'p';
meshFileNameP = 'Mesh Files\PlateDiffusion.msh';
MeshP = BuildMesh_GMSH(meshFileNameP, fieldP, nsd, config_dir, progress_on);
% type of material per element
MeshP.MatList = zeros(MeshP.ne, 1, 'int8');
% assign material type to elements
MeshP.MatList(:) = 1;

%%%% porosity field
if contains(Control.PMmodel, 'UPN')
    fieldN = 'n';
    meshFileNameN = 'Mesh Files\PlateDiffusion.msh';
    MeshN = BuildMesh_GMSH(meshFileNameN, fieldN, nsd, config_dir, progress_on);
    % type of material per element
    MeshN.MatList = zeros(MeshN.ne, 1, 'int8');
    % assign material type to elements
    MeshN.MatList(:) = 1;
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
T = 1000;
BC.fixed_p = MeshP.right_nodes; 
% fixed DOF values
BC.fixed_p_value = @(t) T*ones(length(BC.fixed_p),1);
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
BC.fluxNodes = MeshP.left_nodes;
% prescribed flux [m/s]
BC.flux = -5;
Flux = BC.flux/(length(BC.fluxNodes) - 1);
BC.fluxValue = Flux*ones(size(BC.fluxNodes));

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
Control.tend = 10;   % final simulation time

% Newmark method
Control.beta = 0.7;
Control.gamma = 0.7;
Control.theta = 0.7;

%% Plot data
% DOF to plot graphs
Control.plotu = 2;
Control.plotp = 2;

end