function [Material, MeshU, MeshP, MeshN, BC, Control] = Plate_Dynamic_Komijani(config_dir, progress_on)
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

%% Material properties - Komijani (2019)
% elasticity modulus [GPa]
Material.E = 14.516e-3;
% Poisson's ratio
Material.nu = 0.3;
% porous media permeability [m2/GPa s]
Material.kf = 1.0194e3;
% dynamic viscosity [GPa s]
Material.mu = 1e-12;
% intrinsic permeability [m2]
Material.k = Material.kf * Material.mu;
% fluid bulk modulus [GPa]
Material.Kf = 2.1;
% solid bulk modulus [GPa]
Material.Ks = 1e11;
% material porosity
Material.n = 0.3;
% Biot's coefficient
Material.alpha = 1;
% fluid density [10^9 kg/m3]
Material.rho_f = 1000e-9;
% solid density [10^9 kg/m3]
Material.rho_s = 2000e-9;
% average density of the medium
Material.rho = Material.n*Material.rho_f + (1-Material.n)*Material.rho_s;
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
Material.constLaw = 'PlaneStrain';

%% Mesh parameters
if progress_on
    disp([num2str(toc),': Building Mesh...']);
end

% mesh type
% 'Manual': 1D mesh
% 'Gmsh': 2D mesh, input file from GMSH
MeshType = 'Gmsh';

switch MeshType
    case 'Manual'
        % number of space dimensions
        nsd = 1;
        % number of elements
        ne = 10;
        % column size [m]
        L = 10;
        %%%% solid displacement field
        typeU = 'L3';
        MeshU = Build1DMesh(nsd, ne, L, typeU);
        %%%% fluid pressure field
        typeP = 'L2';
        MeshP = Build1DMesh(nsd, ne, L, typeP);
        %%%% porosity field
        if contains(Control.PMmodel, 'UPN')
            typeN = 'L2';
            MeshN = Build1DMesh(nsd, ne, L, typeN);
        else
            MeshN = [];
        end
    case 'Gmsh'
        % Version 2 ASCII
        % number of space dimensions
        nsd = 2;
        %%%% displacement field
        fieldU = 'u';
        meshFileNameU = 'Mesh Files\Plate_3x3_Q9.msh';
        MeshU = BuildMesh_GMSH(meshFileNameU, fieldU, nsd, config_dir, progress_on);
        %%%% pressure field
        fieldP = 'p';
        meshFileNameP = 'Mesh Files\Plate_3x3_Q4.msh';
        MeshP = BuildMesh_GMSH(meshFileNameP, fieldP, nsd, config_dir, progress_on);
        %%%% porosity field
        if contains(Control.PMmodel, 'UPN')
            fieldN = 'n';
            meshFileNameN = 'Mesh Files\Plate_3x3_Q4.msh';
            MeshN = BuildMesh_GMSH(meshFileNameN, fieldN, nsd, config_dir, progress_on);
        else
            MeshN = [];
        end
end

%% Dirichlet BCs - solid
% displacement u=0 at bottom (y) and right (x)
BC.fixed_u = [MeshU.right_dofx; MeshU.bottom_dofy];
% fixed DOF values
BC.fixed_u_value = zeros(length(BC.fixed_u),1);
% free nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid
% pressure p=0 at boundaries
BC.fixed_p = [MeshP.top_dof; MeshP.bottom_dof; MeshP.left_dof; MeshP.right_dof];
% fixed DOF values
BC.fixed_p_value = zeros(length(BC.fixed_p),1);
% free nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% traction interpolation (needed for traction applied in wells); 1 - true, 0 - false
BC.tractionInterp = 0;
% prescribed traction [GN/m2]
BC.tractionTop = -1e-5;
BC.tractionLeft = -5e-5;

lefttopnode = MeshU.left_nodes(MeshU.coords(MeshU.left_nodes,2) == max(MeshU.coords(:,2)));
index_left = MeshU.left_nodes ~= lefttopnode;
index_top   = MeshU.top_nodes   ~= lefttopnode;
BC.tractionNodes = [MeshU.left_nodes(index_left);  MeshU.top_nodes(index_top); lefttopnode];

ForceTop = BC.tractionTop * max(MeshU.coords(:,1))/((length(MeshU.top_nodes) - 1)/2);
ForceLeft = BC.tractionLeft * max(MeshU.coords(:,2))/((length(MeshU.left_nodes) - 1)/2);

BC.tractionForce = [ForceLeft*ones(size(MeshU.left_nodes(index_left))), zeros(size(MeshU.left_nodes(index_left))); % left side nodes
    zeros(size(MeshU.top_nodes(index_top))), ForceTop*ones(size(MeshU.top_nodes(index_top))); % top side nodes
    ForceLeft*1/2, ForceTop*1/2]; % top right node

% Q9 elements for displacement field
for n = 1:length(BC.tractionForce)
    if any(BC.tractionNodes(n) == MeshU.top_nodes)
        if any(BC.tractionNodes(n) == MeshU.conn(:,1:4),'all')  % then node is a corner node
            BC.tractionForce(n,2) = ForceTop/3;
        else % then node is a midside node
            BC.tractionForce(n,2) = ForceTop*2/3;
        end
    end
end

% Q9 elements for displacement field
for n = 1:length(BC.tractionForce)
    if any(BC.tractionNodes(n) == MeshU.left_nodes)
        if any(BC.tractionNodes(n) == MeshU.conn(:,1:4),'all')  % then node is a corner node
            BC.tractionForce(n,1) = ForceLeft/3;
        else % then node is a midside node
            BC.tractionForce(n,1) = ForceLeft*2/3;
        end
    end
end

% find the nodes in the top left and right corners
righttopnode  = MeshU.right_nodes(MeshU.coords(MeshU.right_nodes,2) == max(MeshU.coords(:,2)));
leftbottomnode = MeshU.left_nodes(MeshU.coords(MeshU.left_nodes,2) == min(MeshU.coords(:,2)));

lefttopnodeindex = find(BC.tractionNodes == lefttopnode);
righttopnodeindex = find(BC.tractionNodes == righttopnode);
leftbottomnodeindex = find(BC.tractionNodes == leftbottomnode);

BC.tractionForce(lefttopnodeindex,2) = BC.tractionForce(lefttopnodeindex,2)/2;
BC.tractionForce(righttopnodeindex,2) = BC.tractionForce(righttopnodeindex,2)/2;

BC.tractionForce(lefttopnodeindex,1) = BC.tractionForce(lefttopnodeindex,1)/2;
BC.tractionForce(leftbottomnodeindex,1) = BC.tractionForce(leftbottomnodeindex,1)/2;

% point loads [GN]
BC.pointLoad = [];

% body force [GN/m3]
BC.b = @(x)[];  

%% Neumann BCs - fluid
% distributed flux [m3/s]
BC.pointFlux = [];

% point flux [m/s]
BC.fluxNodes = [];

% flux source [m3/s/m3]
BC.s = @(x)[]; 

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
Control.dt = 1e-4;  % time step
Control.tend = 0.1;   % final simulation time

% Newmark method
Control.beta = 0.7;
Control.gamma = 0.7;
Control.theta = 0.7;

%% Plot graphs
% DOF to plot graphs
Control.plotu = 121*2; % dof y of node 121 (x = 1.5m, y = 1.5m)
Control.plotp = 81; % dof of node 81 (x = 1.5m, y = 1.5m)

% Fine mesh
% Control.plotu = 681*2; % dof y of node 121 (x = 1.5m, y = 1.5m)
% Control.plotp = 361; % dof of node 81 (x = 1.5m, y = 1.5m)

end