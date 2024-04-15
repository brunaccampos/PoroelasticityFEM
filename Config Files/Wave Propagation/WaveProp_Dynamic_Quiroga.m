function [Material, MeshU, MeshP, MeshN, BC, Control] = WaveProp_Dynamic_Quiroga(~, progress_on,~,~)
% Wave propagation in 2D
% Configuration File
% ------------------------------------------------------------------------
% Based on Zienkiewicz (1982) model for dynamic case
% ------------------------------------------------------------------------
% Assumptions/conventions:
% - stress is positive for tension
% - boundary condition for force is based on total stress
% ------------------------------------------------------------------------

%% Poroelasticity model
% Options:  Tr1_Biot_UP -------- Biot model (u-p), transient
%           Tr2_Spanos_UPN ----- Spanos model (u-p-n), transient
%           Tr3_Spanos_UP ------ Spanos model (u-p), dynamic, implicit
%                                   porosity perturbation equation
%           Dyn1_Biot_UP -------- Biot model (u-p), dynamic
%           Dyn2_Spanos_UPN ----- Spanos model (u-p-n), dynamic
%           Dyn3_Spanos_UP ------ Spanos model (u-p), dynamic, implicit
%                                   porosity perturbation equation
%           Dyn4_Biot_UPU ------- Biot model (u-p-U), dynamic
%           Dyn5_Spanos_UPU ----- Spanos model (u-p-U), dynamic, implicit
%                                   porosity perturbation equation
Control.PMmodel = 'Dyn5_Spanos_UPU';

%% Material properties - Quiroga-Goode (2005)
% Poisson's ratio
Material.nu = 0.2;
% dynamic viscosity [Pa s]
Material.mu = 1e-3;
% intrinsic permeability [m2]
Material.k = 1e-13;
% porous media permeability [m2/Pa s]
Material.kf = Material.k / Material.mu;
% Biot's coefficient
Material.alpha = 0.78;
% fluid bulk modulus [Pa]
Material.Kf = 2.2e9;
% solid bulk modulus [Pa]
Material.Ks = 33e9;
% fluid bulk viscosity [Pa s]
Material.xif = 2.8e-3;
% material porosity
Material.eta0 = 0.25;
% shear modulus [Pa]
Material.G = 4.9e9;
% elasticity modulus [Pa]
Material.E = 2 * Material.G * (1 + Material.nu);
% 1/Q (related to storage coefficient)
Material.Minv = (Material.alpha - Material.eta0)/Material.Ks + Material.eta0/Material.Kf;
% fluid density [kg/m3]
Material.rho_f = 1000;
% solid density [kg/m3]
Material.rho_s = 2650;
% average density of the medium
Material.rho = Material.eta0*Material.rho_f + (1-Material.eta0)*Material.rho_s;
% added mass [kg/m3]
Material.rho12 = -83;

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
% n = Material.Ks/Material.Kf; % upper limit

% modified storage coefficient (Muller, 2019)
Mstarinv = Material.Minv - (1-n)*(Material.alpha - Material.eta0)/Material.Ks; 
Mstar = 1/Mstarinv;

Material.deltaF = (Material.alpha - Material.eta0) * Material.eta0 * Mstar * n / Material.Ks;
Material.deltaS = (Material.alpha - Material.eta0) * Material.eta0 * Mstar /Material.Kf;

%% Mesh parameters
if progress_on
    disp([num2str(toc),': Building Mesh...']);
end

% location of initial node [m] [x0;y0;z0]
coord0 = [0;0;0];
% number of space dimensions
nsd = 2;
% size of domain [m] [Lx;Ly;Lz]
L = [15; 15];
% number of elements in each direction [nex; ney; nez]
ne = [100; 100];

%%%% displacement mesh
% element type ('Q4')
typeU = 'Q4';
% variable field ('u', 'p', 'n')
fieldU = 'u';
MeshU = BuildMesh_structured(nsd, coord0, L, ne, typeU, fieldU, progress_on);

%%%% pressure mesh
% element type ('Q4')
typeP = 'Q4';
% variable field ('u', 'p', 'n')
fieldP = 'p';
MeshP = BuildMesh_structured(nsd, coord0, L, ne, typeP, fieldP, progress_on);

%%%% porosity mesh
if contains(Control.PMmodel, 'UPN')
    % element type ('Q4')
    typeN = 'Q4';
    % variable field ('u', 'p', 'n')
    fieldN = 'n';
    MeshN = BuildMesh_structured(nsd, coord0, L, ne, typeN, fieldN, progress_on);
else
    MeshN = [];
end
        
%% Dirichlet BCs - solid
% central node
node = find(MeshU.coords(:,1) == 7.5 & MeshU.coords(:,2) == 7.5);
% central node y DOF
BC.fixed_u = node*2;
% period [s]
t0 = 1e-3;
% fixed DOF values
BC.fixed_u_value = @(t) (-t0/(2*pi)*cos(2*pi*(t)/t0) + t0/(8*pi)*cos(4*pi*(t)/t0) + 3*t0/8/pi).*(t<t0);
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
BC.pointLoad = [];

% distributed load [GN/m2]
BC.tractionNodes = [];

% body force [GN/m3]
BC.bs = @(x,t)[];  
BC.bf = @(x,t)[];

%% Neumann BCs - fluid
% point flux [m/s]
BC.pointFlux = [];

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
% 1 = uncoupled problem (elasticity, heat transfer, etc)
% 0 = coupled problem (Biot, Spanos model)
Control.uncoupled = 0; 

% plot analytical solution (valid for 1D problems with Material.Minv == 0)
Control.plotansol = 0; % 1 = true; 0 = false

%% Time step controls
Control.dt = 1e-5;  % time step
Control.tend = 3e-3;   % final simulation time

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

Control.depthplot = 7.5; % fixed coordinate
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