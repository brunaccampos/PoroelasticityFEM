function [Material, MeshU, MeshP, MeshN, BC, Control] = SoilSlope(config_dir, progress_on,~,~)
% 2D simulation of soil slope (Wu, 2019)
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
Control.PMmodel = 'Dyn_BT_UPU';

%% Material properties - Wu (2019)
% elasticity modulus [GPa]
Material.M(1).E = 60e-3;
% Poisson's ratio
Material.M(1).nu = 0.4;
% hydraulic permeability [m/s]
Material.M(1).kh = 4e-10;
% gravitational acceleration [m/s2]
Material.M(1).g = 9.81;
% fluid density [Gkg/m3]
Material.M(1).rhof = 1000e-9;
% porous media permeability [m2/Pa s]
Material.M(1).kf = Material.M(1).kh/Material.M(1).rhof/Material.M(1).g;
% dynamic viscosity [GPa s]
Material.M(1).muf = 1e-12;
% intrinsic permeability [m2]
Material.M(1).k = Material.M(1).kf * Material.M(1).muf;
% fluid bulk modulus [GPa]
Material.M(1).Kf = 2;
% solid bulk modulus [GPa]
Material.M(1).Ks = 50;
% material porosity
Material.M(1).eta0 = 0.2;
% Biot's coefficient
Material.M(1).alpha = 1;
% solid density [Gkg/m3]
Material.M(1).rhos = 2000e-9;
% average density of the medium
Material.M(1).rho = Material.M(1).eta0*Material.M(1).rhof + (1-Material.M(1).eta0)*Material.M(1).rhos;
% 1/Q (related to storage coefficient)
Material.M(1).Minv = (Material.M(1).alpha - Material.M(1).eta0)/Material.M(1).Ks + Material.M(1).eta0/Material.M(1).Kf;
% fluid bulk viscosity [GPa s]
Material.M(1).xif = 2.8e-12; % (Quiroga-Goode, 2005)

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

% GMSH file Version 2 ASCII
% number of space dimensions
nsd = 2;
%%%% displacement field
fieldU = 'u';
meshFileNameU = 'Mesh Files\SoilSlopeT3finer.msh';
MeshU = BuildMesh_GMSH(meshFileNameU, fieldU, nsd, config_dir, progress_on);
% type of material per element
MeshU.MatList = zeros(MeshU.ne, 1, 'int8');
% assign material type to elements
MeshU.MatList(:) = 1;

%%%% pressure field
fieldP = 'p';
meshFileNameP = 'Mesh Files\SoilSlopeT3finer.msh';
MeshP = BuildMesh_GMSH(meshFileNameP, fieldP, nsd, config_dir, progress_on);
% type of material per element
MeshP.MatList = zeros(MeshP.ne, 1, 'int8');
% assign material type to elements
MeshP.MatList(:) = 1;

%%%% porosity field
MeshN = [];

%% Load strip 
% select nodes in the interval [(12,14) (16,14)] where the load is applied
coord1 = 12;
coord2 = 16;
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
% displacement u=0 at bottom (x,y), left (x), and right (x)
BC.fixed_u = [MeshU.left_dofx; MeshU.right_dofx; MeshU.bottom_dofx; MeshU.bottom_dofy];
% fixed DOF values
BC.fixed_u_value = @(t) zeros(length(BC.fixed_u),1);
% free nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid displacement
% displacement u=0 at bottom (x,y), left (x), and right (x)
BC.fixed_uf = [MeshU.left_dofx; MeshU.right_dofx; MeshU.bottom_dofx; MeshU.bottom_dofy];
% fixed DOF values
BC.fixed_uf_value = @(t) zeros(length(BC.fixed_uf),1);
% free nodes
BC.free_uf = setdiff(MeshU.DOF, BC.fixed_uf);

%% Dirichlet BCs - fluid
% boundaries for p=0
% nodesP_set1 = [4, 5, 59:75]';
% nodesP_set2 = [3, 4, 42:58]';
% nodesP_set1 = [4, 5, 113:147]'; % fine
% nodesP_set2 = [3, 4, 78:112]'; % fine
nodesP_set1 = [4, 5, 221:291]'; % finer
nodesP_set2 = [3, 4, 150:220]'; % finer
% pressure p=0 at top
BC.fixed_p = unique([MeshP.top_nodes; nodesP_set1; nodesP_set2]);
% fixed DOF values
BC.fixed_p_value = @(t) zeros(length(BC.fixed_p),1);
% free nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% prescribed traction [GN/m2]
BC.traction = -3e-3;
BC.tractionNodes = nodesU;
Force = BC.traction * 1/((length(BC.tractionNodes) - 1)/2);
BC.tractionForce = zeros(length(BC.tractionNodes),2);

% T6 elements for displacement field
if contains(MeshU.type, 'T6')
    for n = 1:length(BC.tractionForce)
        if any(BC.tractionNodes(n) == MeshU.conn(:,1:3),'all') % then node is a corner node
            BC.tractionForce(n,:) = [0, Force/3];
        else % then node is a midside node
            BC.tractionForce(n,:) = [0, Force*2/3];
        end
    end
else
    BC.tractionForce = Force*[zeros(size(BC.tractionNodes)), ones(size(BC.tractionNodes))/2];
end

% find the nodes in the left and right corners of the strip where load is
% applied
lefttopnode = find(MeshU.coords(BC.tractionNodes,1) == 12);
righttopnode  = find(MeshU.coords(BC.tractionNodes,1) == 16);

BC.tractionForce(lefttopnode,2) = BC.tractionForce(lefttopnode,2)/2;
BC.tractionForce(righttopnode,2) = BC.tractionForce(righttopnode,2)/2;

% point loads [GN]
BC.pointLoad = @(t)[];

% body force [GN/m3]
BC.b = @(x,t)[];  

%% Neumann BCs - fluid
% distributed flux [m/s]
% impervious at bottom, left, and right
BC.fluxNodes = [MeshP.left_dof; MeshP.right_dof; MeshP.bottom_dof];
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
Control.nqU = 3;
Control.nqP = 3;

%% Frequency domain
Control.freqDomain = 0;  % 1 = true; 0 = false

%% Analytical solution
% plot analytical solution (valid for 1D problems)
Control.plotansol = 0; % 1 = true; 0 = false

%% Time step controls
Control.dt = 1e-4;  % time step
Control.tend = 1;   % final simulation time

% Newmark method
Control.beta = 0.7;
Control.gamma = 0.7;
Control.theta = 0.7;
Control.lambda = 0.7;

%% Plot data
% DOF to plot graphs
Control.plotu = 5*2; 
Control.plotp = 5; 

% Plot in a row (all nodes at x = 10m)
Control.fixedDepthPlotON = 1; % 0: false, 1: true

Control.depthplot = 14; % fixed coordinate
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