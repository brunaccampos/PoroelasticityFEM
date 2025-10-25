% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function [Material, MeshU, MeshP, MeshN, BC, Control] = ManufacturedSolutionL3(~, progress_on, ~, nelements)
% Manufactured solution for L3 element mesh size convergence study
% ux = x^5 - x^4

%% Poroelasticity model
% use default for running in main function
Control.PMmodel = 'Tr_BT_UP';

%% Material properties
% elasticity modulus [Pa]
Material.M(1).E = 2230;
% Poisson's ratio
Material.M(1).nu = 0.3;
% porous media permeability [m2/Pa s]
Material.M(1).kf = 0;
% 1/Q (related to storage coefficient)
Material.M(1).Minv = 0;
% Biot's coefficient
Material.M(1).alpha = 0;

% thickness 
% 1D: cross sectional area [m2]
% 2D: out of plane thickness [m]
Material.t = 1;

% constititive law - 'PlaneStress' or 'PlaneStrain'
Material.constLaw = 'PlaneStress';

%% Mesh parameters
if progress_on
    disp([num2str(toc),': Building Mesh...']);
end
% number of space dimensions
nsd = 1;
% number of elements
ne = nelements;
% column size [m]
L = 2;
        
% solid displacement field
typeU = 'L3';
fieldU = 'u';
MeshU = Build1DMesh(nsd, ne, L, typeU, fieldU);

% fluid pressure field
typeP = 'L2';
fieldP = 'p';
MeshP = Build1DMesh(nsd, ne, L, typeP, fieldP);

MeshN = [];
        
%% Find nodes for prescribed BCs
% find top and bottom nodes for displacement field
BC.top_node_u = find(MeshU.coords == max(MeshU.coords));
BC.bottom_node_u = find(MeshU.coords == min(MeshU.coords));

% find top and bottom nodes for pressure field
BC.top_node_p = find(MeshP.coords == max(MeshP.coords));
BC.bottom_node_p = find(MeshP.coords == min(MeshP.coords));

%% Dirichlet BCs - solid
BC.ux = @(x) x.^5 - x.^4;
BC.dudx = @(x) 5*x.^4 - 4*x.^3;
% column vector of prescribed displacement dof
BC.fixed_u_dof1 = BC.top_node_u;
BC.fixed_u_dof2 = BC.bottom_node_u;
BC.fixed_u = [BC.fixed_u_dof1; BC.fixed_u_dof2];
% prescribed displacement
BC.fixed_u_value = zeros(length(BC.fixed_u),1);
BC.fixed_u_value = @(t) BC.ux(MeshU.coords(BC.fixed_u));
% free displacement nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid
% prescribed pressure
BC.fixed_p = 1:MeshP.nDOF;
BC.fixed_p_value = @(t) zeros(length(BC.fixed_p),1);
% free pressure nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% column vector of prescribed traction nodes
BC.tractionNodes = [];

% body force
BC.b = @(x,t) - Material.M(1).E * (20 * x.^3 - 12 * x.^2);

% point load [N]
BC.pointLoad = @(t)[];

%% Neumann BCs - fluid
% point flux [m/s]
BC.pointFlux = @(t)[];

% distributed flux [m/s]
BC.fluxNodes = [];

% flux source
BC.s = @(x,t)[]; 

%% Quadrature order
Control.nqU = 3;
Control.nqP = 2;

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

% analytical solution
Control.pan_symb = @(x,t) zeros(MeshP.nDOF,1);
Control.p_an = @(t) Control.pan_symb(MeshP.coords,t);
Control.uan_symb = @(x,t) x.^5 - x.^4;
Control.u_an = @(t) Control.uan_symb(MeshU.coords,t);

%% Time step controls
Control.dt = 1;  % time step
Control.tend = 1; % final simulation time

% Beta method
% beta = 1 Backward Euler; beta = 0.5 Crank-Nicolson
Control.beta = 1; 

%% Plot data
% DOF to plot graphs
Control.plotu = 1;
Control.plotp = 1;

end
