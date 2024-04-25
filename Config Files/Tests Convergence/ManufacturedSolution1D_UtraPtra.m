function [Material, MeshU, MeshP, MeshN, BC, Control] = ManufacturedSolution1D_UtraPtra(~, progress_on,~,~)
% ------------------------------------------------------------------------
% Manufactured solution for L3 element mesh size convergence study
% u = sin(xt)
% p = sin(xt)
% ------------------------------------------------------------------------
% Adapted from https://github.com/GCMLab (Acknowledgements: Bruce Gee)
% ------------------------------------------------------------------------

%% Poroelasticity
% Biot's coefficient
Material.alpha = 0;
% 1/Q (related to storage coefficient)
Material.Minv = 0;
% poroelasticity model
Control.PMmodel = 'Tr1_Biot_UP';

% thickness 
% 1D: cross sectional area [m2]
% 2D: out of plane thickness [m]
Material.t = 1;

% constititive law - 'PlaneStress' or 'PlaneStrain'
% Note: use 'PlaneStrain' for 1D or 2D poroelasticity
Material.constLaw = 'PlaneStress';

%% Material properties
% elasticity modulus [Pa]
Material.E = 1;
% Poisson's ratio
Material.nu = 0.3;
% porous media permeability [m2/Pa s]
Material.kf = 1;

%% Mesh parameters
if progress_on
    disp([num2str(toc),': Building Mesh...']);
end

% location of initial node [m] [x0;y0;z0]
coord0 = [0;0;0];
% number of space dimensions
nsd = 1;
% size of domain [m] [Lx;Ly;Lz]
L = 1;
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
BC.ux = @(x,t) sin(x*t);
% column vector of prescribed displacement dof
BC.fixed_u_dof1 = MeshU.left_nodes;
BC.fixed_u_dof2 = MeshU.right_nodes;
BC.fixed_u = [BC.fixed_u_dof1; BC.fixed_u_dof2];
% prescribed displacement
BC.fixed_u_value = zeros(length(BC.fixed_u),1);
BC.fixed_u_value = @(t) BC.ux(MeshU.coords(BC.fixed_u),t);
% free displacement nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid
BC.p = @(x,t) sin(x*t);
% column vector of prescribed displacement dof
BC.fixed_p_dof1 = MeshP.left_nodes;
BC.fixed_p_dof2 = MeshP.right_nodes;
BC.fixed_p = [BC.fixed_p_dof1; BC.fixed_p_dof2];
% prescribed pressure
BC.fixed_p_value = zeros(length(BC.fixed_p),1);
BC.fixed_p_value = @(t) BC.p(MeshP.coords(BC.fixed_p),t);
% free pressure nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% column vector of prescribed traction nodes
BC.tractionNodes = [];

% body force
BC.b = @(x,t) Material.E * t^2 * sin(x*t);

% point load [N]
BC.pointLoad = @(t)[];

%% Neumann BCs - fluid
% point flux [m/s]
BC.pointFlux = @(t)[];

% distributed flux [m/s]
BC.fluxNodes = [];

% flux source
BC.s = @(x,t) Material.kf * t^2 * sin(x*t); 

%% Quadrature order
Control.nqU = 2;
Control.nqP = 2;

%% Frequency domain
Control.freqDomain = 0;  % 1 = true; 0 = false

%% Analytical solution
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
Control.pan_symb = @(x,t) sin(x*t);
Control.p_an = @(t) Control.pan_symb(MeshP.coords,t);

%% Time step controls
Control.dt = 1e-3;  % time step
Control.tend = 1;

Control.beta = 1; % beta-method time discretization -- beta = 1 Backward Euler; beta = 0.5 Crank-Nicolson

%% Plot data
Control.plotu = round(MeshU.nn/2);
Control.plotp = round(MeshP.nn/2);

end