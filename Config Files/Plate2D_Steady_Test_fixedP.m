function [Material, MeshU, MeshP, MeshN, BC, Control] = Plate2D_Steady_Test_fixedP(config_dir, progress_on)
% Column Consolidation 2D simulation
% Configuration File
% Based on Korsawe (2006) model
% ------------------------------------------------------------------------
% Assumptions/conventions:
% - stress is positive for tension
% - boundary condition for force is based on total stress
% - only solid acceleration is considered (undrained condition; no motions
% of the fluid relative to the solid skeleton can occur)
% - solid grains and fluid are incompressible
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

%% Material properties - Komijani (2019)
% shear modulus [GPa]
Material.G = 2540;
% Poisson's ratio
Material.nu = 0.3;
% elasticity modulus [GPa]
Material.E = 2 * Material.G * (1 + Material.nu);
% porous media permeability [m2/GPa s]
Material.kf = 0;
% dynamic viscosity [GPa s]
Material.mu = 1e-12;
% intrinsic permeability [m2]
Material.k = Material.kf * Material.mu;
% Biot's coefficient
Material.alpha = 0;
% fluid bulk modulus [GPa]
Material.Kf = 3000e-3;
% solid bulk modulus [GPa]
Material.Ks = 36000e-3;
% fluid bulk viscosity [GPa s]
Material.xif = 2.8e-12; % (Quiroga-Goode, 2005)
% material porosity
Material.n = 0.19;
% 1/Q (related to storage coefficient)
% Material.Minv = (Material.alpha - Material.n)/Material.Ks + Material.n/Material.Kf;
Material.Minv = 0;

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

% mesh type
% 'Manual': 1D mesh
% 'Gmsh': 2D mesh, input file from GMSH
MeshType = 'Manual';

switch MeshType
    case 'Manual'
    % Manual 2D mesh
        MeshU.nsd = 2; % number of spatial directions
        MeshU.nn = 9; % number of nodes
        MeshU.ne = 4; % number of elements
        MeshU.type = 'Q4'; % element type
        MeshU.field = 'u'; % field type
        MeshU.nne = 4; % nodes per element
        MeshU.nDOFe = MeshU.nne*MeshU.nsd; % DOFs per element
        MeshU.nDOF = MeshU.nn*MeshU.nsd; % total number of DOFs
        for sd = 1:MeshU.nsd
            MeshU.DOF(:,sd) = (sd : MeshU.nsd : (MeshU.nDOF-(MeshU.nsd-sd)))';
        end
        MeshU.coords = zeros(MeshU.nn, 2); % nodal coordinates
        % x coordinates
        MeshU.coords(:,1) = [0; 1; 1; 0; 0.5; 1; 0.5; 0; 0.5];
        % y coordinates
        MeshU.coords(:,2) = [0; 0; 1; 1; 0; 0.5; 1; 0.5; 0.5];
        % elements connectivity
        MeshU.conn = [1, 5, 9, 8;
                    5, 2, 6, 9;
                    8, 9, 7, 4;
                    9, 6, 3 ,7]; 
                
        % mesh for fluid field
        MeshP.nsd = 2; % number of spatial directions
        MeshP.nn = 9; % number of nodes
        MeshP.ne = 4; % number of elements
        MeshP.type = 'Q4'; % element type
        MeshP.field = 'p'; % field type
        MeshP.nne = 4; % nodes per element
        MeshP.nDOFe = MeshP.nne; % DOFs per element
        MeshP.nDOF = MeshP.nn; % total number of DOFs
        MeshP.DOF = (1:MeshP.nDOF).'; % DOFs
        MeshP.coords = zeros(MeshP.nn, 2); % nodal coordinates
        % x coordinates
        MeshP.coords(:,1) = [0; 1; 1; 0; 0.5; 1; 0.5; 0; 0.5];
        % y coordinates
        MeshP.coords(:,2) = [0; 0; 1; 1; 0; 0.5; 1; 0.5; 0.5];
        % elements connectivity
        MeshP.conn = [1, 5, 9, 8;
                    5, 2, 6, 9;
                    8, 9, 7, 4;
                    9, 6, 3 ,7]; 

        % mesh for porosity field
        MeshN = [];
% ------------------------------------------------------------------------

    case 'Gmsh'
        % Version 2 ASCII
        % number of space dimensions
        nsd = 2;
        %%%% displacement field
        fieldU = 'u';
        meshFileNameU = 'Mesh Files\PlateDiffusion.msh';
        MeshU = BuildMesh_GMSH(meshFileNameU, fieldU, nsd, config_dir, progress_on);
        %%%% pressure field
        fieldP = 'p';
        meshFileNameP = 'Mesh Files\PlateDiffusion.msh';
        MeshP = BuildMesh_GMSH(meshFileNameP, fieldP, nsd, config_dir, progress_on);
        %%%% porosity field
        if contains(Control.PMmodel, 'UPN')
            fieldN = 'n';
            meshFileNameN = 'Mesh Files\PlateDiffusion.msh';
            MeshN = BuildMesh_GMSH(meshFileNameN, fieldN, nsd, config_dir, progress_on);
        else
            MeshN = [];
        end
end

%% Initial conditions
% displacement
BC.initU = [];

% pressure
BC.initP = [];

%% Dirichlet BCs - solid
% ---------- test 1: all boundary displacements fixed in x (represent rigid body
% motion)
BC.fixed_u1 = [1; 3; 5; 7; 9; 11; 13; 15];
BC.fixed_u2 = [2; 4; 6; 8; 10; 12; 14; 16];
BC.fixed_u = [BC.fixed_u1; BC.fixed_u2];

% ---------- test 2: all boundary displacemnts fixed in y (represent rigid body
% motion)
% BC.fixed_u = [MeshU.left_dofy; MeshU.right_dofy; MeshU.bottom_dofy; MeshU.top_dofy];

% ---------- test 3: linear displacement in x (represent constant strain
% state)


% ---------- test 4: linear displacement in y (represent constant strain
% state)


% fixed DOF values
BC.fixed_u_value1 = ones(length(BC.fixed_u1),1);
BC.fixed_u_value2 = zeros(length(BC.fixed_u2),1);
BC.fixed_u_value = [BC.fixed_u_value1; BC.fixed_u_value2];
% free nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid
% pressure fixed at all nodes
BC.fixed_p = (1:MeshP.nDOF);
% fixed DOF values
BC.fixed_p_value = zeros(length(BC.fixed_p),1);
% free nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
BC.tractionNodes = [];

% % prescribed traction [GN/m2]
% BC.traction = 1e-5;
% BC.tractionNodes = MeshU.left_nodes;
% Force = BC.traction * max(MeshU.coords(:,2))/((length(MeshU.left_nodes) - 1)/2);
% BC.tractionForce = zeros(length(BC.tractionNodes),2);
% 
% % Q9 elements for displacement field
% for n = 1:length(BC.tractionForce)
%     if any(BC.tractionNodes(n) == MeshU.conn(:,1:4),'all') % then node is a corner node
%         BC.tractionForce(n,:) = [Force/3, 0];
%     else % then node is a midside node
%         BC.tractionForce(n,:) = [Force*2/3, 0];
%     end
% end
% 
% % find the nodes in the left top and bottom corners
% lefttopnode = find(MeshU.coords(BC.tractionNodes,2) == max(MeshU.coords(:,2)));
% leftbottomnode  = find(MeshU.coords(BC.tractionNodes,2) == min(MeshU.coords(:,2)));
% 
% BC.tractionForce(lefttopnode,1) = BC.tractionForce(lefttopnode,1)/2;
% BC.tractionForce(leftbottomnode,1) = BC.tractionForce(leftbottomnode,1)/2;

% point loads [GN]
BC.pointLoad = [];

% body force [GN/m3]
BC.b = @(x)[];  

%% Neumann BCs - fluid
% distributed flux [m/s]
% impervious at bottom, left, and right
% BC.fluxNodes = [MeshP.left_dof; MeshP.right_dof; MeshP.bottom_dof];
BC.fluxNodes = [];
BC.fluxValue = zeros(length(BC.fluxNodes),1);

% point flux [m3/s]
BC.pointFlux = [];

% flux source [m3/s/m3]
BC.s = @(x)[]; 

%% Quadrature order
Control.nqU = 2;
Control.nqP = 2;

%% Solution parameters
% tag used for computing analytical solution
% 1 = uncoupled problem (elasticity, heat transfer, etc)
% 0 = coupled problem (Biot, Spanos model)
Control.uncoupled = 0; 

% basic time step controls
Control.dt = 1e-2;  % time step
Control.tend = 10;   % final simulation time

Control.beta = 1; % beta-method time discretization -- beta = 1 Backward Euler; beta = 0.5 Crank-Nicolson

% DOF to plot graphs
Control.plotu = 1; % dof y of node 63 (x = 0.05m, y = 3m)
Control.plotp = 1; % dof of node 35 (x = 0.05m, y = 3m)

% plot analytical solution (valid for 1D problems with Material.Minv == 0)
Control.plotansol = 0; % 1 = true; 0 = false

% solve in the frequency domain
Control.freqDomain = 0;  % 1 = true; 0 = false

end