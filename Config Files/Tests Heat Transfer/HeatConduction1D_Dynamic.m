function [Material, MeshU, MeshP, MeshN, BC, Control] = HeatConduction1D_Dynamic(~, progress_on,~,~)
% Heat conduction in 1D
% Configuration File
% ------------------------------------------------------------------------

%% Material properties
% porous media permeability [m2/Pa s]
Material.kf = 1;
% 1/Q (related to storage coefficient)
Material.Minv = 1;
% Poroelasticity model
Control.PMmodel = 'Dyn1_Biot_UP';

% material density [kg/m3]
Material.rho = 0;
% fluid density [kg/m3]
Material.rhof = 0;
% elasticity modulus [Pa]
Material.E = 0;
% Poisson's ratio
Material.nu = 0;
% Biot's coefficient
Material.alpha = 0;
% poroelasticity model
Control.Biotmodel = 1;

% lumped mass matrix - 0: false, 1: true
Material.lumpedMass = 0;

% thickness
% 1D: cross sectional area [m2]
% 2D: out of plane thickness [m]
Material.t = 1;

% constititive law - 'PlaneStress' or 'PlaneStrain'
% Note: use 'PlaneStrain' for 1D or 2D poroelasticity
Material.constLaw = 'PlaneStress';

%% Mesh parameters
if progress_on
    disp([num2str(toc),': Building Mesh...']);
end

% location of initial node [m] [x0;y0;z0]
coord0 = [0;0;0];
% number of space dimensions
nsd = 1;
% size of domain [m] [Lx;Ly;Lz]
L = 20;
% number of elements in each direction [nex; ney; nez]
ne = 10;

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
% column vector of prescribed displacement dof
BC.fixed_u = 1:MeshU.nDOF;
% prescribed displacement for each dof [u1; u2; ...] [m]
BC.fixed_u_value = @(t) zeros(length(BC.fixed_u),1);
% free nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid
% pressure p=0 at the top
BC.fixed_p = MeshP.left_nodes;
% fixed DOF values
BC.fixed_p_value = @(t) zeros(length(BC.fixed_p),1);
% free nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% distributed traction [N/m2]
BC.tractionNodes = [];

% point loads [N]
BC.pointLoad = @(t)[];

% body force [N/m3]
BC.b = @(x,t)[];

%% Neumann BCs - fluid
% point flux [m3/s]
BC.pointFluxValue = -1;
BC.pointFluxNodes = MeshP.right_nodes;
BC.pointFlux = zeros(MeshP.nDOF,1);
BC.pointFlux(BC.pointFluxNodes) = BC.pointFluxValue;

% distributed flux [m/s]
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
Control.nqU = 2;
Control.nqP = 2;

%% Frequency domain
Control.freqDomain = 0;  % 1 = true; 0 = false

%% Analytical solution
% plot analytical solution (valid for 1D problems with Material.Minv == 0)
Control.plotansol = 0; % 1 = true; 0 = false

%% Time step controls
Control.dt = 0.01;  % time step
Control.tend = 1;   % final simulation time

% Newmark method
Control.beta = 0.7;
Control.gamma = 0.7;
Control.theta = 0.7;

%% Plot data
% DOF to plot graphs
Control.plotu = 1;
Control.plotp = 1;

end