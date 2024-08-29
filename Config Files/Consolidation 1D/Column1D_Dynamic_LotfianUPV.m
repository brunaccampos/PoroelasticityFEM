function [Material, MeshU, MeshP, MeshN, BC, Control] = Column1D_Dynamic_LotfianUPV(~, progress_on,~,~)
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
Control.PMmodel = 'Dyn6_Biot_UPV';

%% Material properties - Komijani (2019)
% elasticity modulus [GPa]
Material.M(1).E = 14.516e-3;
% Poisson's ratio
Material.M(1).nu = 0.3;
% porous media permeability [m2/GPa s]
Material.M(1).kf = 1.0194e3;
% dynamic viscosity [GPa s]
Material.M(1).muf = 1e-12;
% intrinsic permeability [m2]
Material.M(1).k = Material.M(1).kf * Material.M(1).muf;
% fluid bulk modulus [GPa]
Material.M(1).Kf = 2.1;
% solid bulk modulus [GPa]
Material.M(1).Ks = 1e11;
% material porosity
Material.M(1).eta0 = 0.3;
% Biot's coefficient
Material.M(1).alpha = 1;
% fluid density [10^9 kg/m3]
Material.M(1).rhof = 1000e-9;
% solid density [10^9 kg/m3]
Material.M(1).rhos = 2000e-9;
% average density of the medium
Material.M(1).rho = Material.M(1).eta0*Material.M(1).rhof + (1-Material.M(1).eta0)*Material.M(1).rhos;
% 1/Q (related to storage coefficient)
Material.M(1).Minv = (Material.M(1).alpha - Material.M(1).eta0)/Material.M(1).Ks + Material.M(1).eta0/Material.M(1).Kf;
% fluid bulk viscosity [GPa s]
Material.M(1).xif = 2.8e-12; % (Quiroga-Goode, 2005)

% thickness
% 1D: cross sectional area [m2]
% 2D: out of plane thickness [m]
Material.t = 0.1;

% constititive law - 'PlaneStress' or 'PlaneStrain'
% Note: use 'PlaneStrain' for 1D or 2D poroelasticity
Material.constLaw = 'PlaneStrain';

%% Spanos material parameters
% porosity effective pressure coefficient (Spanos, 1989)
% n = 0; % lower limit
n = 1; % return to Biot
% n = Material.M(1).Ks/Material.M(1).Kf; % upper limit

% modified storage coefficient (Muller, 2019)
Mstarinv = Material.M(1).Minv - (1-n)*(Material.M(1).alpha - Material.M(1).eta0)/Material.M(1).Ks; 
Mstar = 1/Mstarinv;

Material.M(1).deltaf = (Material.M(1).alpha - Material.M(1).eta0) * Material.M(1).eta0 * Mstar * n / Material.M(1).Ks;
Material.M(1).deltas = (Material.M(1).alpha - Material.M(1).eta0) * Material.M(1).eta0 * Mstar /Material.M(1).Kf;

%% Mesh parameters
if progress_on
    disp([num2str(toc),': Building Mesh...']);
end

% location of initial node [m] [x0;y0;z0]
coord0 = [0;0;0];
% number of space dimensions
nsd = 1;
% size of domain [m] [Lx;Ly;Lz]
L = 10;
% number of elements in each direction [nex; ney; nez]
ne = 400;

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

MeshN = [];

%% Dirichlet BCs - solid
% displacement u=0 at the bottom
BC.fixed_u = MeshU.right_nodes;
BC.fixed_u_value = @(t) zeros(length(BC.fixed_u),1);
% free displacement nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid velocity
% velocity u=0 at the bottom
BC.fixed_ufdot = MeshU.right_nodes;
BC.fixed_ufdot_value = @(t) zeros(length(BC.fixed_ufdot),1);
% free displacement nodes
BC.free_ufdot = setdiff(MeshU.DOF, BC.fixed_ufdot);

%% Dirichlet BCs - fluid pressure
% pressure p=0 at the top
BC.fixed_p = MeshP.left_nodes;
BC.fixed_p_value = @(t) zeros(length(BC.fixed_p),1);
% free pressure nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% point load [GN]
BC.pointLoadValue = 3e-6;
BC.pointLoadNodes = MeshU.left_nodes;
BC.pointLoad = zeros(MeshU.nDOF,1);
BC.pointLoad(BC.pointLoadNodes) = BC.pointLoadValue;

% distributed load [GN/m2]
BC.tractionNodes = [];

% body force [GN/m3]
BC.bs = @(x,t) [];

BC.bf = @(x,t) [];

%% Neumann BCs - fluid
% point flux [m/s]
BC.pointFlux = @(t)[];

% distributed flux [m3/s]
BC.fluxNodes = [];

% flux source [m3/s/m3]
BC.s = @(x,t) [];

%% Quadrature order
Control.nqU = 3;
Control.nqP = 3;

%% Frequency domain
Control.freqDomain = 0;  % 1 = true; 0 = false

%% Analytical solution
% plot analytical solution (valid for 1D problems)
Control.plotansol = 0; % 1 = true; 0 = false

%% Time step controls
Control.dt = 1e-4;  % time step [s]
Control.tend = 1;   % final simulation time [s]

% Newmark method
Control.beta = 0.7;
Control.gamma = 0.7;
Control.theta = 0.7;
Control.lambda = 0.7;

%% Plot data
% DOF to plot graphs
Control.plotu = round(length(MeshU.coords)/2);
Control.plotp = round(length(MeshP.coords)/2);

% Plot in a row
Control.fixedDepthPlotON = 0; % 0: false, 1: true

% Nodes to plot in a row (all nodes for 1D case)
Control.ploturow = MeshU.DOF;
Control.plotprow = MeshP.DOF;

end