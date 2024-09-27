function [Material, MeshU, MeshP, MeshN, BC, Control] = Plate2D_Dynamic_Pulse(~, progress_on,~,~)
% Pulse propagation in 2D simulation
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
Control.PMmodel = 'Dyn_dCS_UPU';

%% Material properties - Berea Sandstone (Detournay, 1993, p.26)
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

% thickness
% 1D: cross sectional area [m2]
% 2D: out of plane thickness [m]
Material.t = 1;

% constititive law - 'PlaneStress' or 'PlaneStrain'
% Note: use 'PlaneStrain' for 1D or 2D poroelasticity
Material.constLaw = 'PlaneStrain';

% lumped mass matrix - 0: false, 1: true
Material.lumpedMass = 0;

% lumped damping matrix - 0: false, 1: true
Material.lumpedDamping = 0;

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
nsd = 2;
% size of domain [m] [Lx;Ly;Lz]
L = [12; 12];
% number of elements in each direction [nex; ney; nez]
ne = [240; 240];

%%%% displacement mesh
% element type ('Q4')
typeU = 'Q4';
% variable field ('u', 'p', 'n')
fieldU = 'u';
MeshU = BuildMesh_structured(nsd, coord0, L, ne, typeU, fieldU, progress_on);
% type of material per element
MeshU.MatList = zeros(MeshU.ne, 1, 'int8');
% assign material type to elements
MeshU.MatList(:) = 1;

%%%% pressure mesh
% element type ('Q4')
typeP = 'Q4';
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
    typeN = 'Q4';
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
% central node
node = find(MeshU.coords(:,1) == 6 & MeshU.coords(:,2) == 6);
% central node y DOF
BC.fixed_u = node*2;
% amplitude [GPa]
P0 = 100e-6;
% frequency [Hz]
f = 1e3;
BC.fixed_u_value = @(t) (P0*(sin(pi*(t)*f) - 0.5*sin(2*pi*(t)*f))).*(t<1/f);
% BC.fixed_u_value = @(t) (P0*(sin(2*pi*(t)*f) - 0.5*sin(4*pi*(t)*f))).*(t<1/f);
% free displacement nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid
BC.fixed_p = [MeshP.right_dof; MeshP.left_dof; MeshP.top_dof; MeshP.bottom_dof];
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
% plot analytical solution (valid for 1D problems with Material.Minv == 0)
Control.plotansol = 0; % 1 = true; 0 = false

%% Time step controls
Control.dt = 1e-5;  % time step
Control.tend = 2e-3;   % final simulation time

% Newmark method
Control.beta = 0.7;
Control.gamma = 0.7;
Control.theta = 0.7;
Control.lambda = 0.7;

% HHT method (-1/3 < alpha < 0)
Control.alpha = 0;

% adaptive time step (optional)
% Control.dtmin = 1e-3; % minimum time step
% Control.tlim = 1; % limit to use dtmin

% ramp load option (optional); uses tlim from adaptive time step
% NOTE: only declare if true
% Control.rampLoad = 1;

%% Plot data
% DOF to plot graphs
nodeU = find(MeshU.coords(:,1) == 6 & MeshU.coords(:,2) == 5);
nodeP = find(MeshP.coords(:,1) == 6 & MeshP.coords(:,2) == 5);
Control.plotu = nodeU*2; % dof y at x = 6m, y = 5m
Control.plotp = nodeP; % dof y x = 6m, y = 5m

% Plot synthetics
Control.plotSyntheticsON = 1; % 0: false, 1: true

% Plot in a row
Control.fixedDepthPlotON = 1; % 0: false, 1: true

Control.depthplot = 6; % fixed coordinate
Control.depthDir = 2; % 1 = fixed y, vary x --- 2 = fixed x, vary y

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