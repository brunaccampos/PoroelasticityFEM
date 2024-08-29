function [Material, MeshU, MeshP, MeshN, BC, Control] = Column1D_Dynamic_Gajo(~, progress_on,~,~)
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
% Tr1_Biot_UP -------- Biot model (u-p), transient, u = solid
%                      displacement, p = fluid pressure
% Tr2_Spanos_UPN ----- Spanos model (u-p-n), transient, u = solid
%                      displacement, p = fluid pressure, n = porosity
% Tr3_Spanos_UP ------ Spanos model (u-p), dynamic, u = solid displacement,
%                      p = fluid pressure
% Dyn1_Biot_UP ------- Biot model (u-p), dynamic, u = solid displacement, 
%                      p = fluid pressure
% Dyn2_Spanos_UPN ---- Spanos model (u-p-n), dynamic, u = solid displacement,
%                      p = fluid pressure, n = porosity
% Dyn3_Spanos_UP ----- Spanos model (u-p), dynamic, u = solid displacement, 
%                      p = fluid pressure
% Dyn4_Biot_UPU ------ Biot model (u-p-U), dynamic, u = solid displacement,
%                      p = fluid pressure, U = fluid displacement
% Dyn5_Spanos_UPU ---- Spanos model (u-p-U), dynamic, u = solid displacement,
%                      p = fluid pressure, U = fluid displacement
% Dyn6_Biot_UPV ------ Biot model (u-p-v), dynamic, u = solid displacement,
%                      p = fluid pressure, v = fluid velocity
% Dyn7_Biot_UPW ------ Biot model (u-p-w), dynamic, u = solid displacement,
%                      p = fluid pressure, w = relative fluid velocity
Control.PMmodel = 'Dyn4_Biot_UPU';

%% Material properties - Tasiopoulou (2015)
% elasticity modulus [GPa]
Material.M(1).E = 0.2321e2;
% Poisson's ratio
Material.M(1).nu = 0.171;
% porous media permeability [m2/GPa s]
Material.M(1).kf = 0.1480e4;
% dynamic viscosity [GPa s]
Material.M(1).muf = 1e-12;
% intrinsic permeability [m2]
Material.M(1).k = Material.M(1).kf * Material.M(1).muf;
% fluid bulk modulus [GPa]
Material.M(1).Kf = 2.2;
% solid bulk modulus [GPa]
Material.M(1).Ks = 36;
% material porosity
Material.M(1).eta0 = 0.18;
% Biot's coefficient
Material.M(1).alpha = 0.677;
% fluid density [10^9 kg/m3]
Material.M(1).rhof = 1000e-9;
% solid density [10^9 kg/m3]
Material.M(1).rhos = 2660e-9;
% average density of the medium
Material.M(1).rho = Material.M(1).eta0*Material.M(1).rhof + (1-Material.M(1).eta0)*Material.M(1).rhos;
% 1/Q (related to storage coefficient)
Material.M(1).Minv = (Material.M(1).alpha - Material.M(1).eta0)/Material.M(1).Ks + Material.M(1).eta0/Material.M(1).Kf;
% fluid bulk viscosity [GPa s]
Material.M(1).xif = 2.8e-12; % (Quiroga-Goode, 2005)

% thickness 
% 1D: cross sectional area [m2]
% 2D: out of plane thickness [m]
Material.t = 1;

% constititive law - 'PlaneStress' or 'PlaneStrain'
% Note: use 'PlaneStrain' for 1D or 2D poroelasticity
Material.constLaw = 'PlaneStrain';

% lumped mass matrix - 0: false, 1: true
Material.lumpedMass = 1;

% lumped damping matrix - 0: false, 1: true
Material.lumpedDamping = 1;

%% Spanos material parameters
% porosity effective pressure coefficient (Spanos, 1989)
% n = 0; % lower limit
n = 1; % return to Biot
% n = Material.M(1).Ks/Material.M(1).Kf; % upper limit

% modified storage coefficient (Muller, 2019)
Mstarinv = Material.M(1).Minv - (1-n)*(Material.M(1).alpha - Material.M(1).eta0)/Material.M(1).Ks; 
Mstar = 1/Mstarinv;

% porosity equation coefficients
Material.M(1).deltaf = (Material.M(1).alpha - Material.M(1).eta0) * Material.M(1).eta0 * Mstar * n / Material.M(1).Ks;
Material.M(1).deltas = (Material.M(1).alpha - Material.M(1).eta0) * Material.M(1).eta0 * Mstar / Material.M(1).Kf;

%% Mesh parameters
if progress_on
    disp([num2str(toc),': Building Mesh...']);
end

% location of initial node [m] [x0;y0;z0]
coord0 = [0;0;0];
% number of space dimensions
nsd = 1;
% size of domain [m] [Lx;Ly;Lz]
L = 0.5;
% number of elements in each direction [nex; ney; nez]
ne = 100;

%%%% displacement mesh
% element type ('Q4')
typeU = 'L3';
% variable field ('u', 'p', 'n')
fieldU = 'u';
MeshU = BuildMesh_structured(nsd, coord0, L, ne, typeU, fieldU, progress_on);
% type of material per element
MeshU.MatList = zeros(MeshU.ne, 1, 'int8');
% assign material type to elements
MeshU.MatList(:) = 1;

%%%% pressure mesh
% element type ('Q4')
typeP = 'L2';
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

%% Dirichlet BCs - solid
% displacement u=0 at the bottom
BC.fixed_u1 = MeshU.right_nodes;
% u = t at the top
BC.fixed_u2 = MeshU.left_nodes;
% fixed DOFs
BC.fixed_u = [BC.fixed_u1; BC.fixed_u2];
% fixed DOF values
BC.fixed_u_value = @(t) [zeros(length(BC.fixed_u1),1); ones(length(BC.fixed_u2),1)];
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
BC.pointLoad = @(t)[];

% distributed load [GN/m2]
BC.tractionNodes = [];

% body force [GN/m3]
BC.b = @(x,t)[];  

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
% plot analytical solution (valid for 1D problems)
Control.plotansol = 0; % 1 = true; 0 = false

%% Time step controls
Control.dt = 1.4e-6;  % time step [s]
Control.tend = 1e-3;   % final simulation time [s]

% Newmark method
Control.beta = 0.6;
Control.gamma = 0.6;
Control.theta = 0.6;
Control.lambda = 0.6;

%% Plot data
% DOF to plot graphs
Control.plotu = find(MeshU.coords == 0.1); % x = 0.1m
Control.plotp = find(MeshP.coords == 0.1); % x = 0.1m

end