function [Material, MeshU, MeshP, MeshN, BC, Control] = ManufacturedSolution1D_UstePtra(~, progress_on,~,~)
% ------------------------------------------------------------------------
% Manufactured solution for element mesh size convergence study
% u = sin(x) - steady
% p = BVP - transient
% ------------------------------------------------------------------------
% Adapted from https://github.com/GCMLab (Acknowledgements: Bruce Gee)
% ------------------------------------------------------------------------

%% Poroelasticity
% Biot's coefficient
Material.M(1).alpha = 0;
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
% porous media permeability [m2/Pa s]
Material.M(1).kf = 0.002;
% 1/Q (related to storage coefficient)
Material.M(1).Minv = 1;
% elasticity modulus [Pa]
Material.M(1).E = 1;
% Poisson's ratio
Material.M(1).nu = 0.3;

%% Mesh parameters
if progress_on
    disp([num2str(toc),': Building Mesh...']);
end

% location of initial node [m] [x0;y0;z0]
coord0 = [0;0;0];
% number of space dimensions
nsd = 1;
% size of domain [m] [Lx;Ly;Lz]
L = 2;
% number of elements in each direction [nex; ney; nez]
ne = 128;

%%%% displacement mesh
% element type ('Q4')
typeU = 'L4';
% variable field ('u', 'p', 'n')
fieldU = 'u';
MeshU = BuildMesh_structured(nsd, coord0, L, ne, typeU, fieldU, progress_on);
% type of material per element
MeshU.MatList = zeros(MeshU.ne, 1, 'int8');
% assign material type to elements
MeshU.MatList(:) = 1;

%%%% pressure mesh
% element type ('Q4')
typeP = 'L3';
% variable field ('u', 'p', 'n')
fieldP = 'p';
MeshP = BuildMesh_structured(nsd, coord0, L, ne, typeP, fieldP, progress_on);
% type of material per element
MeshP.MatList = zeros(MeshP.ne, 1, 'int8');
% assign material type to elements
MeshP.MatList(:) = 1;

%%%% porosity mesh
if contains(Control.PMmodel, 'UPN')
    % element type ('Q4')
    typeN = 'L2';
    % variable field ('u', 'p', 'n')
    fieldN = 'n';
    MeshN = BuildMesh_structured(nsd, coord0, L, ne, typeN, fieldN, progress_on);
    % type of material per element
    MeshN.MatList = zeros(MeshN.ne, 1, 'int8');
    % assign material type to elements
    MeshN.MatList(:) = 1;
else
    MeshN = [];
end

%% Initial conditions
% pressure
BC.initP = 50*ones(MeshP.nDOF,1);

%% Dirichlet BCs - solid
BC.ux = @(x,t) sin(x);
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
% pressure p=0 at the top
BC.fixed_p = [MeshP.left_nodes; MeshP.right_nodes];
% fixed DOF values
BC.fixed_p_value = @(t) [100; 0];
% free nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% column vector of prescribed traction nodes
BC.tractionNodes = [];

% body force
BC.b = @(x,t) sin(x);

% point load [N]
BC.pointLoad = @(t)[];

%% Neumann BCs - fluid
% point flux [m/s]
BC.pointFlux = @(t)[];

% distributed flux [m/s]
BC.fluxNodes = [];

% flux source
BC.s = @(x,t) [];

%% Quadrature order
Control.nqU = 3;
Control.nqP = 3;

%% Frequency domain
Control.freqDomain = 0;  % 1 = true; 0 = false

%% Time step controls
Control.dt = 1e-2;  % time step
Control.tend = 1;

Control.beta = 1; % beta-method time discretization -- beta = 1 Backward Euler; beta = 0.5 Crank-Nicolson

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
Control.uan_symb = @(x,t) sin(x);
Control.u_an = @(t) Control.uan_symb(MeshU.coords,t);
aux=0;
syms x
N=1000;
for k=1:N
    aux = aux + (1/k)*exp(-Material.kf*k^2*pi()^2*Control.tend)*sin(k*pi()*x);
end
T0vec = BC.fixed_p_value(0);
steady = (T0vec(2) - T0vec(1))*x/max(MeshP.coords) + T0vec(1);

% solution in p
Control.pan_symb = @(t) steady - (T0vec(1)/pi()) * aux; % transient solution
aux=0;
x = MeshP.coords;
N=1000;
for k=1:N
    aux = aux + (1/k)*exp(-Material.kf*k^2*pi()^2*Control.tend)*sin(k*pi()*x);
end
steady = (T0vec(2) - T0vec(1))*x/max(MeshP.coords) + T0vec(1);
Control.p_an = @(t) steady - (T0vec(1)/pi()) * aux; % transient solution

%% Plot data
Control.plotu = round(MeshU.nn/2);
Control.plotp = round(MeshP.nn/2);

end