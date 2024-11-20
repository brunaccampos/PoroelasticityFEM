function [Material, MeshU, MeshP, MeshN, BC, Control] = WaveProp_Dynamic_Tian(~, progress_on,~,~)
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
Control.PMmodel = 'Dyn_dCS_UPU';

%% Material properties - Tian (2023)
vps = 3000; %  P wave solid velocity [m/s]
vss = 1732; % S wave solid velocity [m/s]
vpf = 1500; % P wave fluid velocity [m/s]
eta0 = 0.15; % porosity [-]
rhos = 2588; % solid density [kg/m3]
rhof = 952; % fluid density [kg/m3]
rho12 = 0.5*(eta0-1)*rhof; % added mass [kg/m3]

Ks = rhos*(vps^2-4*vss^2/3); % solid bulk modulus [Pa]
Kf = rhof*vpf^2; % fluid bulk modulus [Pa]
mus = rhos*vss^2; % solid shear modulus [Pa]

e = mus*(9*Ks+8*mus)/6/(Ks+2*mus);
Kd = (4*mus*Ks*(1-eta0))/(4*mus+3*eta0*Ks);
mud = e*mus*(1-eta0)/(e+eta0*mus);
Kk = Kd + (1-Kd/Ks)^2/(eta0/Kf+(1-eta0/Ks-Kd/Ks^2));
L = 1/((1-eta0-Kk/Ks)/Ks+eta0/Kf);
K = (1-Kk/Ks)*L;
H = (1-Kk/Ks)^2*L+Kk+4*mud/3;

A = H-2*K*eta0+L*eta0^2-2*mud;
N = mud;
Q = K*eta0-L*eta0^2;
R = L*eta0^2;

alpha = (Q/R+1)*eta0;

% Poisson's ratio
Material.M(1).nu = 0.2;
% dynamic viscosity [Pa s]
Material.M(1).muf = 1e-3;
% intrinsic permeability [m2]
Material.M(1).k = 1e-13;
% porous media permeability [m2/Pa s]
Material.M(1).kf = Material.M(1).k / Material.M(1).muf;
% Biot's coefficient
Material.M(1).alpha = alpha;
% fluid bulk modulus [Pa]
Material.M(1).Kf = Kf;
% solid bulk modulus [Pa]
Material.M(1).Ks = Ks;
% fluid bulk viscosity [Pa s]
Material.M(1).xif = 2.8e-3;
% material porosity
Material.M(1).eta0 = eta0;
% shear modulus [Pa]
Material.M(1).mu = N;
% elasticity modulus [Pa]
Material.M(1).E = 2 * Material.M(1).mu * (1 + Material.M(1).nu);
% 1/Q (related to storage coefficient)
Material.M(1).Minv = (Material.M(1).alpha - Material.M(1).eta0)/Material.M(1).Ks + Material.M(1).eta0/Material.M(1).Kf;
% fluid density [kg/m3]
Material.M(1).rhof = rhof;
% solid density [kg/m3]
Material.M(1).rhos = rhos;
% average density of the medium
Material.M(1).rho = Material.M(1).eta0*Material.M(1).rhof + (1-Material.M(1).eta0)*Material.M(1).rhos;
% added mass [kg/m3]
Material.M(1).rho12 = rho12;

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
nsd = 2;
% size of domain [m] [Lx;Ly;Lz]
L = [5000; 5000];
% number of elements in each direction [nex; ney; nez]
ne = [500; 500];

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
node = find(MeshU.coords(:,1) == 2500 & MeshU.coords(:,2) == 2500);
% central node y DOF
BC.fixed_u = node*2;
% peak frequency [Hz]
f = 10;
% peak location [s]
t0 = 1/f;
% fixed DOF values
BC.fixed_u_value = @(t) (1-2*(pi*f*(t-t0)).^2) .* exp(-(pi*f*(t-t0)).^2);
% t = 0:0.001:0.5;
% plot(t, BC.fixed_u_value(t));
% free displacement nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid displacement
% displacement prescribed on the left and right
BC.fixed_uf = [];
BC.fixed_uf_value = @(t) zeros(length(BC.fixed_uf),1);
% free displacement nodes
BC.free_uf = setdiff(MeshU.DOF, BC.fixed_uf);

%% Dirichlet BCs - fluid
BC.fixed_p = [];
% fixed DOF values
BC.fixed_p_value = @(t) zeros(length(BC.fixed_p),1);
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
Control.plotu = node*2; % dof y of node 242 (x = 7.5m, y = 7.5m)
Control.plotp = node; % dof y of node 177 (x = 7.5m, y = 7.5m)

% Plot synthetics
Control.plotSyntheticsON = 1; % 0: false, 1: true

% Plot in a row
Control.fixedDepthPlotON = 1; % 0: false, 1: true

Control.depthplot = 2600; % fixed coordinate
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

% number of parallel solvers
Control.parallel = 8;

end