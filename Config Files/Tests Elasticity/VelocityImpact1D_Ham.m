function [Material, MeshU, MeshP, MeshN, BC, Control] = VelocityImpact1D_Ham(~, progress_on,~,~)
% Velocity impact for 1D elasticity
% Configuration File
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
Control.PMmodel = 'Dyn1_Biot_UP';

%% Material properties - Idesman (2011)
% elasticity modulus [Pa]
Material.E = 200e9;
% average density of the medium
Material.rho = 8000;
% Poisson's ratio
Material.nu = 0.0;

% porous media permeability [m2/GPa s]
Material.kf = 0;
% Biot's coefficient
Material.alpha = 0;
% 1/Q (related to storage coefficient)
Material.Minv = 0;
% fluid density [10^9 kg/m3]
Material.rho_f = 0;

% lumped mass matrix - 0: false, 1: true
Material.lumpedMass = 1;

% thickness 
% 1D: cross sectional area [m2]
% 2D: out of plane thickness [m]
Material.t = 0.01;

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
L = 0.5;
% number of elements in each direction [nex; ney; nez]
ne = 100;

%%%% displacement mesh
% element type ('Q4')
typeU = 'L2';
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
% displacement u=0 at right
BC.fixed_u1 = MeshU.right_nodes;
% displacement u=t at left (v=1)
BC.fixed_u2 = MeshU.left_nodes;
% fixed DOFs
BC.fixed_u = [BC.fixed_u1; BC.fixed_u2];
% fixed DOF values
BC.fixed_u_value = @(t) [zeros(length(BC.fixed_u1),1); t*ones(length(BC.fixed_u2),1)];
% free nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid
% pressure p=0 at all boundaries
BC.fixed_p = 1:MeshP.nDOF;
% fixed DOF values
BC.fixed_p_value = @(t) zeros(length(BC.fixed_p),1);
% free nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% prescribed traction [GN/m2]
BC.tractionNodes = [];

% point loads [GN]
BC.pointLoad = [];

% body force [GN/m3]
BC.b = @(x,t)[];  

%% Neumann BCs - fluid
% distributed flux [m3/s]
BC.fluxNodes = [];

% point flux [m/s]
BC.pointFlux = [];

% flux source [m3/s/m3]
BC.s = @(x,t)[]; 

%% Quadrature order
Control.nqU = 3;
Control.nqP = 3;

%% Frequency domain
Control.freqDomain = 0;  % 1 = true; 0 = false

%% Analytical solution
% 1 = uncoupled problem (elasticity, heat transfer, etc)
% 0 = coupled problem (Biot, Spanos model)
Control.uncoupled = 1; 

% plot analytical solution (valid for 1D problems with Material.Minv == 0)
Control.plotansol = 0; % 1 = true; 0 = false

%% Time step controls
Control.dt = 2.5e-8;  % time step
Control.tend = 5e-5;   % final simulation time

% Newmark method
Control.beta = 0.25;
Control.gamma = 0.5;
Control.theta = 0.5;

%% Plot data
% DOF to plot graphs
Control.plotu = find(MeshU.coords == max(MeshU.coords)/2); % middle node
Control.plotp = find(MeshU.coords == max(MeshU.coords)/2); % middle node

end