% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function [Material, MeshU, MeshP, MeshN, BC, Control] = WavePropP_MatLayer2000x2500(config_dir, progress_on,~,~)
% Wave propagation in 2D
% Configuration File
% ------------------------------------------------------------------------
% Based on Zienkiewicz (1982) model for dynamic case
% ------------------------------------------------------------------------
% Assumptions/conventions:
% - stress is positive for tension
% - boundary condition for force is based on total stress
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
Control.PMmodel = 'Dyn_BT_UPW';

%% Material properties - berea sandstone (Detournay, 1993)
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

%% Material properties - shale
mus = 18.4; Ks = 35; eta0 = 0.08; alpha = 0.47; mu = (1-eta0)*mus; K = (1-alpha)*Ks;
% elasticity modulus [GPa]
Material.M(2).E = 9*K*mu/(3*K+mu);
% Poisson's ratio
Material.M(2).nu = (3*K-2*mu)/(2*(3*K+mu));
% intrinsic permeability [m2]
Material.M(2).k = 1e-21;
% dynamic viscosity [GPa s]
Material.M(2).muf = 1e-12;
% porous media permeability [m2/GPa s]
Material.M(2).kf = Material.M(2).k/Material.M(2).muf;
% Biot's coefficient
Material.M(2).alpha = 0.47;
% fluid bulk modulus [GPa]
Material.M(2).Kf = 3.3;
% solid bulk modulus [GPa]
Material.M(2).Ks = 35;
% material porosity
Material.M(2).eta0 = 0.08;
% 1/Q (related to storage coefficient)
Material.M(2).Minv = (Material.M(2).alpha - Material.M(2).eta0)/Material.M(2).Ks + Material.M(2).eta0/Material.M(2).Kf;
% fluid bulk viscosity [GPa s]
Material.M(2).xif = 2.8e-12;
% fluid density [10^9 kg/m3]
Material.M(2).rhof = 1000e-9;
% solid density [10^9 kg/m3]
Material.M(2).rhos = 2600e-9;
% average density of the medium
Material.M(2).rho = Material.M(2).eta0*Material.M(2).rhof + (1-Material.M(2).eta0)*Material.M(2).rhos;

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
% n = Material.M(1).Ks/Material.M(1).Kf; % upper limit

% modified storage coefficient (Muller, 2019) - material 1
Mstarinv = Material.M(1).Minv - (1-n)*(Material.M(1).alpha - Material.M(1).eta0)/Material.M(1).Ks; 
Mstar = 1/Mstarinv;

Material.M(1).deltaf = (Material.M(1).alpha - Material.M(1).eta0) * Material.M(1).eta0 * Mstar * n / Material.M(1).Ks;
Material.M(1).deltas = (Material.M(1).alpha - Material.M(1).eta0) * Material.M(1).eta0 * Mstar /Material.M(1).Kf;

% modified storage coefficient (Muller, 2019) - material 2
Mstarinv = Material.M(2).Minv - (1-n)*(Material.M(2).alpha - Material.M(2).eta0)/Material.M(2).Ks; 
Mstar = 1/Mstarinv;

Material.M(2).deltaf = (Material.M(2).alpha - Material.M(2).eta0) * Material.M(2).eta0 * Mstar * n / Material.M(2).Ks;
Material.M(2).deltas = (Material.M(2).alpha - Material.M(2).eta0) * Material.M(2).eta0 * Mstar /Material.M(2).Kf;

%% Mesh parameters
if progress_on
    disp([num2str(toc),': Building Mesh...']);
end

% Version 2 ASCII
% number of space dimensions
nsd = 2;
%%%% displacement field
fieldU = 'u';
meshFileNameU = 'Mesh Files\LayeredDomain_2000x2500.msh';
MeshU = BuildMesh_GMSH(meshFileNameU, fieldU, nsd, config_dir, progress_on);
% type of material per element
MeshU.MatList = zeros(MeshU.ne, 1, 'int8');
% find nodes correspondent to each material
nodes_mat1 = find(MeshU.coords(:,2) > 900 & MeshU.coords(:,2) < 1100);
nodes_mat2 = setdiff((1:MeshU.nn)', nodes_mat1);
% find elements correspondent to selected nodes
econn_mat1 = ismember(MeshU.conn, nodes_mat1);
elements_mat1 = find(any(econn_mat1,2));
elements_mat2 = setdiff((1:MeshU.ne)', elements_mat1);
% assign material type to elements
MeshU.MatList(elements_mat1) = 1;
MeshU.MatList(elements_mat2) = 2;
% assign material type to nodes
MeshU.MatNodes = zeros(MeshU.ne, 1, 'int8');
MeshU.MatNodes(nodes_mat1) = 1;
MeshU.MatNodes(nodes_mat2) = 2;

%%%% pressure field
fieldP = 'p';
meshFileNameP = 'Mesh Files\LayeredDomain_2000x2500.msh';
MeshP = BuildMesh_GMSH(meshFileNameP, fieldP, nsd, config_dir, progress_on);
% type of material per element
MeshP.MatList = zeros(MeshP.ne, 1, 'int8');
% assign material type to elements
MeshP.MatList(elements_mat1) = 1;
MeshP.MatList(elements_mat2) = 2;

%%%% porosity mesh
MeshN = [];
        
%% Dirichlet BCs - solid
% central node
BC.fixed_u = [MeshU.left_dofy; MeshU.right_dofy; MeshU.top_dofx; MeshU.bottom_dofx];
BC.fixed_u_value = @(t) zeros(length(BC.fixed_u),1);
% free displacement nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid displacement
% displacement prescribed on the left and right
BC.fixed_w = [MeshU.left_dofy];
BC.fixed_w_value = @(t) zeros(length(BC.fixed_w),1);
% free displacement nodes
BC.free_w = setdiff(MeshU.DOF, BC.fixed_w);

%% Dirichlet BCs - fluid
% boundary node
node = find(MeshP.coords(:,1) == 0 & MeshP.coords(:,2) == 1000);
BC.fixed_p = [node; MeshP.top_nodes; MeshP.bottom_nodes; MeshP.right_nodes];
% peak frequency [Hz]
f = 10;
% period [s]
t0 = 1/f;
BC.fixed_p_value = @(t) [(sin(2*pi*(t)*f) - 0.5*sin(4*pi*(t)*f)).*(t<t0); zeros(length(BC.fixed_p)-1,1)];
% free pressure nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% point load [GN]
BC.pointLoad = @(t)[];

% distributed load [GN/m2]
BC.tractionNodes = [];

% body force [GN/m3]
BC.bs = @(x,t)[];  
BC.bf = @(x,t)[];

%% Neumann BCs - fluid
% point flux [m/s]
BC.pointFlux = @(t)[];

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
% plot analytical solution (valid for 1D problems with Material.Minv == 0)
Control.plotansol = 0; % 1 = true; 0 = false

%% Time step controls
Control.dt = 8e-4;  % time step
Control.tend = 8e-1;   % final simulation time

% Newmark method
Control.beta = 0.7;
Control.gamma = 0.7;
Control.theta = 0.7;
Control.lambda = 0.7;

%% Plot data
% DOF to plot graphs
Control.plotu = node*2; % dof y of node where source is applied
Control.plotp = node; % dof of node where source is applied

% Plot synthetics
Control.plotSyntheticsON = 1; % 0: false, 1: true

% Plot in a row
Control.fixedDepthPlotON = 1; % 0: false, 1: true

Control.depthplot = 1000; % fixed coordinate
Control.depthDir = 1; % 1 = fixed y, vary x --- 2 = fixed x, vary y

% node numbering
switch Control.depthDir
    case 1
        rowofnodes_u = find(MeshU.coords(:,2) == Control.depthplot);
        rowofnodes_p = find(MeshP.coords(:,2) == Control.depthplot); 
    case 2
        rowofnodes_u = find(MeshU.coords(:,1) == Control.depthplot); 
        rowofnodes_p = find(MeshP.coords(:,1) == Control.depthplot); 
end

nodes_u = [MeshU.coords(rowofnodes_u,Control.depthDir), rowofnodes_u]; % matrix with node numbering and variable coord
nodes_u_sorted = sortrows(nodes_u); % order in terms of variable coord
Control.ploturow = [nodes_u_sorted(:,2) .* 2 - 1; nodes_u_sorted(:,2) .* 2];

nodes_p = [MeshP.coords(rowofnodes_p,Control.depthDir), rowofnodes_p]; % matrix with node numbering and variable coord
nodes_p_sorted = sortrows(nodes_p); % order in terms of variable coord
Control.plotprow = nodes_p_sorted(:,2);

end