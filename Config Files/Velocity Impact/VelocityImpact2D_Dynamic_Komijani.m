function [Material, MeshU, MeshP, MeshN, BC, Control] = VelocityImpact2D_Dynamic_Komijani(config_dir, progress_on,~,~)
% Column Consolidation 2D simulation
% Configuration File
% Based on Zienkiewicz (1982) model
%
% Assumptions/conventions:
% - stress is positive for tension
% - boundary condition for force is based on total stress
% - only solid acceleration is considered (undrained condition; no motions
% of the fluid relative to the solid skeleton can occur)
% - solid grains and fluid are incompressible

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
Control.PMmodel = 'Dyn1_Biot_UP';

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

% lumped mass matrix - 0: false, 1: true
Material.lumpedMass = 0;

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

% porosity equation coefficients
Material.M(1).deltaf = (Material.M(1).alpha - Material.M(1).eta0) * Material.M(1).eta0 * Mstar * n / Material.M(1).Ks;
Material.M(1).deltas = (Material.M(1).alpha - Material.M(1).eta0) * Material.M(1).eta0 * Mstar / Material.M(1).Kf;

%% Mesh parameters
if progress_on
    disp([num2str(toc),': Building Mesh...']);
end

% Version 2 ASCII
% number of space dimensions
nsd = 2;
%%%% displacement field
fieldU = 'u';
meshFileNameU = 'Mesh Files\VelocityImpactQ9.msh';
MeshU = BuildMesh_GMSH(meshFileNameU, fieldU, nsd, config_dir, progress_on);
% type of material per element
MeshU.MatList = zeros(MeshU.ne, 1, 'int8');
% assign material type to elements
MeshU.MatList(:) = 1;

%%%% pressure field
fieldP = 'p';
meshFileNameP = 'Mesh Files\VelocityImpactQ4.msh';
MeshP = BuildMesh_GMSH(meshFileNameP, fieldP, nsd, config_dir, progress_on);
% type of material per element
MeshP.MatList = zeros(MeshP.ne, 1, 'int8');
% assign material type to elements
MeshP.MatList(:) = 1;

%%%% porosity field
if contains(Control.PMmodel, 'UPN')
    fieldN = 'n';
    meshFileNameN = 'Mesh Files\VelocityImpactQ4.msh';
    MeshN = BuildMesh_GMSH(meshFileNameN, fieldN, nsd, config_dir, progress_on);
    % type of material per element
    MeshN.MatList = zeros(MeshN.ne, 1, 'int8');
    % assign material type to elements
    MeshN.MatList(:) = 1;
else
    MeshN = [];
end

%% Dirichlet BCs - solid
% displacement u=0 at bottom (y), top (y), and right (x)
BC.fixed_u1 = [MeshU.right_dofx; MeshU.top_dofy; MeshU.bottom_dofy];
% u = t at left (x)
BC.fixed_u2 = MeshU.left_dofx;
% fixed DOFs
BC.fixed_u = [BC.fixed_u1; BC.fixed_u2];
% fixed DOF values
BC.fixed_u_value = @(t) [zeros(length(BC.fixed_u1),1); t*ones(length(BC.fixed_u2),1)];
% free nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid
% pressure p=0 at all boundaries
BC.fixed_p = [MeshP.top_dof; MeshP.bottom_dof; MeshP.right_dof; MeshP.left_dof];
% fixed DOF values
BC.fixed_p_value = @(t) zeros(length(BC.fixed_p),1);
% free nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% prescribed traction [GN/m2]
BC.tractionNodes = [];

% point loads [GN]
BC.pointLoad = @(t)[];

% body force [GN/m3]
BC.b = @(x,t)[];  

%% Neumann BCs - fluid
% distributed flux [m3/s]
BC.fluxNodes = [];

% point flux [m/s]
BC.pointFlux = @(t)[];

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
Control.dt = 7e-4;  % time step
Control.tend = 0.3;   % final simulation time

% Newmark method
Control.beta = 0.7;
Control.gamma = 0.7;
Control.theta = 0.7;
Control.lambda = 0.7;

%% Plot data
% DOF to plot graphs
Control.plotu = 98*2-1; % dof x of node 98 (x = 3m, y = 0.05m)
Control.plotp = 54; % dof of node 54 (x = 3m, y = 0.05m)

% Plot synthetics
Control.plotSyntheticsON = 0; % 0: false, 1: true

% Plot in a row
Control.fixedDepthPlotON = 1; % 0: false, 1: true

Control.depthplot = 0.05; % fixed coordinate
tol = 1e-12;
Control.depthDir = 1; % 1 = fixed y, vary x --- 2 = fixed x, vary y

% node numbering
switch Control.depthDir
    case 1
        rowofnodes_u = find(abs(MeshU.coords(:,2) - Control.depthplot) < tol);
        rowofnodes_p = find(abs(MeshP.coords(:,2) - Control.depthplot) < tol); 
    case 2
        rowofnodes_u = find(abs(MeshU.coords(:,1) - Control.depthplot) < tol); 
        rowofnodes_p = find(abs(MeshP.coords(:,1) - Control.depthplot) < tol); 
end

nodes_u = [MeshU.coords(rowofnodes_u,Control.depthDir), rowofnodes_u]; % matrix with node numbering and variable coord
nodes_u_sorted = sortrows(nodes_u); % order in terms of variable coord
Control.ploturow = [nodes_u_sorted(:,2) .* 2 - 1; nodes_u_sorted(:,2) .* 2];

nodes_p = [MeshP.coords(rowofnodes_p,Control.depthDir), rowofnodes_p]; % matrix with node numbering and variable coord
nodes_p_sorted = sortrows(nodes_p); % order in terms of variable coord
Control.plotprow = nodes_p_sorted(:,2);

end