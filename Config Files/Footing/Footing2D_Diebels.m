function [Material, MeshU, MeshP, MeshN, BC, Control] = Footing2D_Diebels(config_dir, progress_on,~,~)
% 2D simulation of footing problem
% Configuration File
% ------------------------------------------------------------------------
% Based on Korsawe (2006) model for transient/quasi-steady case
% ------------------------------------------------------------------------
% Assumptions/conventions:
% - stress is positive for tension
% - boundary condition for force is based on total stress
% - no acceleration terms for solid or fluid
% - solid velocity is neglected
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
Control.PMmodel = 'Tr1_Biot_UP';

%% Material properties - Diebels (1996)
% Lame constant [GPa]
lambda = 5.583e-3;
% shear modulus [GPa]
G = 8.375e-3;
% elasticity modulus [GPa]
Material.M(1).E = G*(3*lambda+2*G)/(lambda+G);
% Poisson's ratio
Material.M(1).nu = 0.2;
% porous media permeability [m2/GPa s]
Material.M(1).kf = 1e3;
% dynamic viscosity [GPa s]
Material.M(1).muf = 1e-12;
% intrinsic permeability [m2]
Material.M(1).k = Material.M(1).kf * Material.M(1).muf;
% material porosity
Material.M(1).eta0 = 0.33;
% Biot's coefficient
Material.M(1).alpha = 1;
% 1/Q (related to storage coefficient)
Material.M(1).Minv = 0;

% % alternative values of Kf, Ks for compressible materials (Boone, 1990)
% % fluid bulk modulus [GPa]
% Material.M(1).Kf = 3;
% % solid bulk modulus [GPa]
% Material.M(1).Ks = 36;
% % 1/Q (related to storage coefficient)
% Material.M(1).Minv = (Material.M(1).alpha - Material.M(1).eta0)/Material.M(1).Ks + Material.M(1).eta0/Material.M(1).Kf;

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
% n = Material.Ks/Material.Kf; % upper limit
 
% porosity equation coefficients
Material.M(1).deltaf = 0;
Material.M(1).deltas = Material.M(1).alpha - Material.M(1).eta0;

% % alternative equations for compressible materials
% % modified storage coefficient (Muller, 2019)
% Mstarinv = Material.M(1).Minv - (1-n)*(Material.M(1).alpha - Material.M(1).eta0)/Material.M(1).Ks; 
% Mstar = 1/Mstarinv;
% 
% % porosity equation coefficients
% Material.M(1).deltaf = (Material.M(1).alpha - Material.M(1).eta0) * Material.M(1).eta0 * Mstar * n / Material.M(1).Ks;
% Material.M(1).deltas = (Material.M(1).alpha - Material.M(1).eta0) * Material.M(1).eta0 * Mstar / Material.M(1).Kf;

%% Mesh parameters
if progress_on
    disp([num2str(toc),': Building Mesh...']);
end

% GMSH file Version 2 ASCII
% number of space dimensions
nsd = 2;
%%%% displacement field
fieldU = 'u';
meshFileNameU = 'Mesh Files\Footing_DiebelsQ9uniformCoarse.msh';
MeshU = BuildMesh_GMSH(meshFileNameU, fieldU, nsd, config_dir, progress_on);
% type of material per element
MeshU.MatList = zeros(MeshU.ne, 1, 'int8');
% assign material type to elements
MeshU.MatList(:) = 1;

%%%% pressure field
fieldP = 'p';
meshFileNameP = 'Mesh Files\Footing_DiebelsQ4uniformCoarse.msh';
MeshP = BuildMesh_GMSH(meshFileNameP, fieldP, nsd, config_dir, progress_on);
% type of material per element
MeshP.MatList = zeros(MeshP.ne, 1, 'int8');
% assign material type to elements
MeshP.MatList(:) = 1;

%%%% porosity field
if contains(Control.PMmodel, 'UPN')
    fieldN = 'n';
    meshFileNameN = 'Mesh Files\Footing_DiebelsQ4uniformCoarse.msh';
    MeshN = BuildMesh_GMSH(meshFileNameN, fieldN, nsd, config_dir, progress_on);
    % type of material per element
    MeshN.MatList = zeros(MeshN.ne, 1, 'int8');
    % assign material type to elements
    MeshN.MatList(:) = 1;
else
    MeshN = [];
end

%% Load strip 
% select nodes in the interval [(0,10) (5,10)] where the load is applied
coord1 = 5;
coord2 = 10;
nu = 1;
for i = 1:length(MeshU.top_nodes)
   if MeshU.coords(MeshU.top_nodes(i),1) >= coord1 &&  MeshU.coords(MeshU.top_nodes(i),1) <= coord2
      nodesU(nu,1) = MeshU.top_nodes(i);
      nu = nu+1;
   end
end

np = 1;
for i = 1:length(MeshP.top_nodes)
   if MeshP.coords(MeshP.top_nodes(i),1) >= coord1 &&  MeshP.coords(MeshP.top_nodes(i),1) <= coord2
      nodesP(np,1) = MeshP.top_nodes(i);
      np = np+1;
   end
end

%% Dirichlet BCs - solid
% displacement u=0 at bottom (y), left (x), and right (x)
BC.fixed_u = [MeshU.left_dofx; MeshU.right_dofx; MeshU.bottom_dofy];
% fixed DOF values
BC.fixed_u_value = @(t) zeros(length(BC.fixed_u),1);
% free nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid
% pressure p=0 at top
BC.fixed_p = setdiff(MeshP.top_nodes, nodesP);
% fixed DOF values
BC.fixed_p_value = @(t) zeros(length(BC.fixed_p),1);
% free nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% prescribed traction [GN/m2]
BC.traction = -15e-6;
BC.tractionNodes = nodesU;
Force = BC.traction * 1/((length(BC.tractionNodes) - 1)/2);
BC.tractionForce = zeros(length(BC.tractionNodes),2);

% Q9 elements for displacement field
for n = 1:length(BC.tractionForce)
    if any(BC.tractionNodes(n) == MeshU.conn(:,1:4),'all') % then node is a corner node
        BC.tractionForce(n,:) = [0, Force/3];
    else % then node is a midside node
        BC.tractionForce(n,:) = [0, Force*2/3];
    end
end

% find the nodes in the left and right corners of the strip where load is
% applied
lefttopnode = find(MeshU.coords(BC.tractionNodes,1) == 5);
righttopnode  = find(MeshU.coords(BC.tractionNodes,1) == max(MeshU.coords(:,1)));

BC.tractionForce(lefttopnode,2) = BC.tractionForce(lefttopnode,2)/2;
BC.tractionForce(righttopnode,2) = BC.tractionForce(righttopnode,2)/2;

% point loads [GN]
BC.pointLoad = @(t)[];

% body force [GN/m3]
BC.b = @(x,t)[];  

%% Neumann BCs - fluid
% distributed flux [m/s]
% impervious at bottom, left, and right
BC.fluxNodes = [MeshP.left_dof; MeshP.right_dof; MeshP.bottom_dof; nodesP];
BC.fluxValue = zeros(length(BC.fluxNodes),1);

% point flux [m3/s]
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
Control.nqU = 2;
Control.nqP = 2;

%% Frequency domain
Control.freqDomain = 0;  % 1 = true; 0 = false

%% Analytical solution
% plot analytical solution (valid for 1D problems)
Control.plotansol = 0; % 1 = true; 0 = false

%% Time step controls
Control.dt = 1e-2;  % time step
Control.tend = 10;   % final simulation time

Control.beta = 1; % beta-method time discretization -- beta = 1 Backward Euler; beta = 0.5 Crank-Nicolson

%% Plot data
% DOF to plot graphs
% Control.plotu = 227*2; % dof y at (x = 5m, y = 10m)
% Control.plotp = 127; % dof at (x = 5m, y = 10m)
% fine mesh
% Control.plotu = 1601*2; % dof y at (x = 5m, y = 5m)
% Control.plotp = 1401; % dof at (x = 5m, y = 5m)
% coarse mesh
Control.plotu = 341*2; % dof y at (x = 5m, y = 5m)
Control.plotp = 261; % dof at (x = 5m, y = 5m)

% Plot in a row (all nodes at x = 10m)
Control.fixedDepthPlotON = 1; % 0: false, 1: true

Control.depthplot = 10; % fixed coordinate
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