function [Material, MeshU, MeshP, MeshN, BC, Control] = Plate2D_KirschTest_BC2(config_dir, progress_on)
% Plate with hole for Kirsch test
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
Control.PMmodel = 'Tr1_Biot_UP';

%% Material properties - Berea Sandstone (Detournay, 1993, p.26)
% elasticity modulus [GPa]
Material.E = 14.4;
% Poisson's ratio
Material.nu = 0.2;
% intrinsic permeability [m2]
Material.k = 1.88e-13;
% dynamic viscosity [GPa s]
Material.mu = 1e-12;
% porous media permeability [m2/GPa s]
Material.kf = Material.k/Material.mu;
% Biot's coefficient
Material.alpha = 0.79;
% fluid bulk modulus [GPa]
Material.Kf = 3.3;
% solid bulk modulus [GPa]
Material.Ks = 36;
% material porosity
Material.n = 0.19;
% 1/Q (related to storage coefficient)
Material.Minv = (Material.alpha - Material.n)/Material.Ks + Material.n/Material.Kf;

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

% GMSH file Version 2 ASCII
% number of space dimensions
nsd = 2;
%%%% displacement field
fieldU = 'u';
meshFileNameU = 'Mesh Files\PlateWithHole_FullGeometryQ4.msh';
MeshU = BuildMesh_GMSH(meshFileNameU, fieldU, nsd, config_dir, progress_on);
%%%% pressure field
fieldP = 'p';
meshFileNameP = 'Mesh Files\PlateWithHole_FullGeometryQ4.msh';
MeshP = BuildMesh_GMSH(meshFileNameP, fieldP, nsd, config_dir, progress_on);
%%%% porosity field
if contains(Control.PMmodel, 'UPN')
    fieldN = 'n';
    meshFileNameN = 'Mesh Files\PlateWithHole_FullGeometryQ4.msh';
    MeshN = BuildMesh_GMSH(meshFileNameN, fieldN, nsd, config_dir, progress_on);
else
    MeshN = [];
end

%% Dirichlet BCs - solid
% displacement
BC.fixed_u = [];
% fixed DOF values
BC.fixed_u_value = @(t) zeros(length(BC.fixed_u),1);
% free nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid
% fixed DOFs
BC.fixed_p = MeshP.DOF;
% fixed DOF values
BC.fixed_p_value = @(t) zeros(length(BC.fixed_p),1);
% free nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% horizontal stress
sigmaH = 10e-3;
% vertical stress
sigmaV = 10e-3;
% corner nodes
topleftnode = MeshU.left_nodes(MeshU.coords(MeshU.left_nodes,2) == max(MeshU.coords(:,2)));
toprightnode = MeshU.right_nodes(MeshU.coords(MeshU.right_nodes,2) == max(MeshU.coords(:,2)));
bottomleftnode = MeshU.left_nodes(MeshU.coords(MeshU.left_nodes,2) == min(MeshU.coords(:,2)));
bottomrightnode = MeshU.right_nodes(MeshU.coords(MeshU.right_nodes,2) == min(MeshU.coords(:,2)));
% index of corner nodes
index_left = MeshU.left_nodes ~= topleftnode & MeshU.left_nodes ~= bottomleftnode;
index_right = MeshU.right_nodes ~= toprightnode & MeshU.right_nodes ~= bottomrightnode;
index_top = MeshU.top_nodes ~= topleftnode & MeshU.top_nodes ~= toprightnode;
index_bottom = MeshU.bottom_nodes ~= bottomleftnode & MeshU.bottom_nodes ~= bottomrightnode;
% prescribed traction nodes
BC.tractionNodes = [MeshU.left_nodes(index_left); MeshU.right_nodes(index_right); ...
    MeshU.top_nodes(index_top); MeshU.bottom_nodes(index_bottom); ...
    topleftnode; toprightnode; bottomleftnode; bottomrightnode];

% prescribed traction value
tractionXvalue = sigmaH*max(MeshU.coords(:,2))/((length(MeshU.right_nodes) - 1));
tractionYvalue = sigmaV*max(MeshU.coords(:,1))/((length(MeshU.top_nodes) - 1));
BC.tractionForce = [tractionXvalue*ones(size(MeshU.left_nodes(index_left))), zeros(size(MeshU.left_nodes(index_left))); % left side nodes
    -tractionXvalue*ones(size(MeshU.right_nodes(index_right))), zeros(size(MeshU.right_nodes(index_right))); % right side nodes
    zeros(size(MeshU.top_nodes(index_top))), -tractionYvalue*ones(size(MeshU.top_nodes(index_top))); % top side nodes
    zeros(size(MeshU.bottom_nodes(index_bottom))), tractionYvalue*ones(size(MeshU.bottom_nodes(index_bottom))); % bottom side nodes
    tractionXvalue*1/2, -tractionYvalue*1/2; % top left node
    -tractionXvalue*1/2, -tractionYvalue*1/2; % top right node
    tractionXvalue*1/2, tractionYvalue*1/2; % bottom left node
    -tractionXvalue*1/2, tractionYvalue*1/2]; % bottom right node

% time dependent
BC.tractionForce = @(t) BC.tractionForce;

% point loads [GN] - initial stress state
BC.pointLoadNodes = 1:MeshU.nDOF;
BC.pointLoad = zeros(MeshU.nDOF, 1);
BC.pointLoad(1:2:end) = sigmaH; % stress in x
BC.pointLoad(2:2:end) = sigmaV; % stress in y

% body force [GN/m3]
BC.b = @(x,t)[];

%% Neumann BCs - fluid
% distributed flux [m3/s]
% impervious at bottom, left, and right
BC.fluxNodes = [];
BC.fluxValue = zeros(length(BC.fluxNodes),1);

% point flux [m/s]
BC.pointFlux = [];

% flux source [m3/s/m3]
BC.s = @(x,t)[];

%% Quadrature order
Control.nqU = 2;
Control.nqP = 2;

%% Frequency domain
Control.freqDomain = 0;  % 1 = true; 0 = false

%% Analytical solution
% 1 = uncoupled problem (elasticity, heat transfer, etc)
% 0 = coupled problem (Biot, Spanos model)
Control.uncoupled = 0;

% plot analytical solution (valid for 1D problems with Material.Minv == 0)
Control.plotansol = 0; % 1 = true; 0 = false

%% Time step controls
Control.dt = 1;  % time step [s]
Control.tend = 1;   % final simulation time [s]

Control.beta = 1; % beta-method time discretization -- beta = 1 Backward Euler; beta = 0.5 Crank-Nicolson

%% Plot data
% DOF to plot graphs
Control.plotu = 9*2; % DOF at well
Control.plotp = 9; % DOF at well

Control.ploturow = MeshU.bottom_dofy;
Control.plotprow = MeshP.bottom_nodes;

end