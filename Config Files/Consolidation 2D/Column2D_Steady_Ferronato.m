% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function [Material, MeshU, MeshP, MeshN, BC, Control] = Column2D_Steady_Ferronato(config_dir, progress_on,~,~)
% Column Consolidation 1D simulation
% Configuration File
% Based on Korsawe (2006) model
% ------------------------------------------------------------------------
% Assumptions/conventions:
% - stress is positive for tension
% - boundary condition for force is based on total stress
% - no acceleration terms for solid or fluid
% - solid velocity is neglected
% - fluid and solid grains are incompressible
% - porosity is constant in space and varies over time
% ------------------------------------------------------------------------
% column top at x=L, column bottom at x=0
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
Control.PMmodel = 'Tr_BT_UP';

%% Material properties - Ferronato (2010)
% shear modulus [GPa]
Material.M(1).mu = 40e-3;
% Lame constant [GPa]
Material.M(1).lambda = 40e-3;
% Poisson's ratio
Material.M(1).nu = Material.M(1).lambda/(2*(Material.M(1).lambda + Material.M(1).mu));
% elasticity modulus [GPa]
Material.M(1).E = 2 * Material.M(1).mu * (1 + Material.M(1).nu);
% gravitational acceleration [m/s2]
Material.M(1).g = 9.81;
% fluid density [10^9 kg/m3]
Material.M(1).rhof = 1000e-9;
% hydraulic conductivity [m/s]
Material.M(1).kh = 1e-5;
% porous media permeability [m2/GPa s]
Material.M(1).kf = Material.M(1).kh/(Material.M(1).rhof * Material.M(1).g);
% dynamic viscosity [GPa s]
Material.M(1).muf = 1e-12;
% intrinsic permeability [m2]
Material.M(1).k = Material.M(1).kf * Material.M(1).muf;
% Biot's coefficient
Material.M(1).alpha = 1;
% fluid bulk modulus [GPa]
Material.M(1).Kf = 1/(4.4e-4) *1e-3;

% fluid bulk viscosity [GPa s]
Material.M(1).xif = 2.8e-12; % (Quiroga-Goode, 2005)
% material porosity
Material.M(1).eta0 = 0.375;
% 1/Q (related to storage coefficient)
Material.M(1).Minv = Material.M(1).eta0/Material.M(1).Kf;

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

% porosity equation coefficients
Material.M(1).deltaf = 0;
Material.M(1).deltas = Material.M(1).alpha - Material.M(1).eta0;

%% Mesh parameters
if progress_on
    disp([num2str(toc),': Building Mesh...']);
end

% Version 2 ASCII
% number of space dimensions
nsd = 2;
%%%% displacement field
fieldU = 'u';
meshFileNameU = 'Mesh Files\Column2DQ9_Ferronato.msh';
MeshU = BuildMesh_GMSH(meshFileNameU, fieldU, nsd, config_dir, progress_on);
% type of material per element
MeshU.MatList = zeros(MeshU.ne, 1, 'int8');
% assign material type to elements
MeshU.MatList(:) = 1;
%%%% pressure field
fieldP = 'p';
meshFileNameP = 'Mesh Files\Column2DQ4_Ferronato.msh';
MeshP = BuildMesh_GMSH(meshFileNameP, fieldP, nsd, config_dir, progress_on);
% type of material per element
MeshP.MatList = zeros(MeshP.ne, 1, 'int8');
% assign material type to elements
MeshP.MatList(:) = 1;
%%%% porosity field
if contains(Control.PMmodel, 'UPN')
    fieldN = 'n';
    meshFileNameN = 'Mesh Files\Column2DQ4_Ferronato.msh';
    MeshN = BuildMesh_GMSH(meshFileNameN, fieldN, nsd, config_dir, progress_on);
    % type of material per element
    MeshN.MatList = zeros(MeshN.ne, 1, 'int8');
    % assign material type to elements
    MeshN.MatList(:) = 1;
else
    MeshN = [];
end

%% Dirichlet BCs - solid
% displacement u=0 at the bottom
BC.fixed_u = [MeshU.left_dofx; MeshU.right_dofx; MeshU.top_dofy];
BC.fixed_u_value = @(t) zeros(length(BC.fixed_u),1);
% free displacement nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid
%   pressure p=0 at the top
BC.fixed_p = (MeshP.bottom_dof);
BC.fixed_p_value = @(t) zeros(length(BC.fixed_p),1);
% free pressure nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% traction interpolation (needed for traction applied in wells); 1 - true, 0 - false
BC.tractionInterp = 0;
% prescribed traction [GN/m2]
BC.traction = 1e-5;
BC.tractionNodes = MeshU.bottom_nodes;
Force = BC.traction * max(MeshU.coords(:,1))/((length(MeshU.top_nodes) - 1)/2);
BC.tractionForce = zeros(length(BC.tractionNodes),2);

% Q9 elements for displacement field
for n = 1:length(BC.tractionForce)
    if any(BC.tractionNodes(n) == MeshU.conn(:,1:4),'all') % then node is a corner node
        BC.tractionForce(n,:) = [0, Force/3];
    else % then node is a midside node
        BC.tractionForce(n,:) = [0, Force*2/3];
    end
end

% find the nodes in the top left and right corners
lefttopnode = find(MeshU.coords(BC.tractionNodes,1) == min(MeshU.coords(:,1)));
righttopnode  = find(MeshU.coords(BC.tractionNodes,1) == max(MeshU.coords(:,1)));

BC.tractionForce(lefttopnode,2) = BC.tractionForce(lefttopnode,2)/2;
BC.tractionForce(righttopnode,2) = BC.tractionForce(righttopnode,2)/2;

% point loads [GN]
BC.pointLoad = @(t)[];

% body force [GN/m3]
BC.b = @(x,t)[];  

%% Neumann BCs - fluid
% point flux [m/s]
BC.fluxNodes = [MeshP.left_dof; MeshP.right_dof; MeshP.top_dof];
BC.fluxValue = zeros(length(BC.fluxNodes),1);

% point flux [m3/s]
BC.pointFlux = @(t)[];

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
% plot analytical solution (valid for 1D problems)
Control.plotansol = 0; % 1 = true; 0 = false

%% Time step controls
Control.dt = 0.1;  % time step [s]
Control.tend = 60;   % final simulation time [s]

Control.beta = 1; % beta-method time discretization -- beta = 1 Backward Euler; beta = 0.5 Crank-Nicolson

%% Plot data
% DOF to plot graphs
Control.plotu = 147*2; % dof y of node 70 (x = 0.033m, y = 7.5m)
Control.plotp = 81; % dof y of node 42 (x = 0.033m, y = 7.5m)

end