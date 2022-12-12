function [Material, MeshU, MeshP, MeshN, BC, Control] = Column2D_Steady_Test_fixedP(config_dir, progress_on)
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
% column top at x=L, column bottom at x=0
% ------------------------------------------------------------------------

%% Poroelasticity model
% 1 - Biot theory
% 0 - Spanos theory (additional porosity equation)
Control.Biotmodel = 1;

%% Material properties - Komijani (2019)
% shear modulus [GPa]
Material.G = 6000e-3;
% Poisson's ratio
Material.nu = 0.2;
% elasticity modulus [GPa]
Material.E = 2 * Material.G * (1 + Material.nu);
% porous media permeability [m2/GPa s]
Material.kf = 2e-2;
% dynamic viscosity [GPa s]
Material.mu = 1e-12;
% intrinsic permeability [m2]
Material.k = Material.kf * Material.mu;
% Biot's coefficient
Material.alpha = 1;
% fluid bulk modulus [GPa]
Material.Kf = 3000e-3;
% solid bulk modulus [GPa]
Material.Ks = 36000e-3;
% fluid bulk viscosity [GPa s]
Material.xif = 2.8e-12; % (Quiroga-Goode, 2005)
% material porosity
Material.n = 0.19;
% 1/Q (related to storage coefficient)
Material.Minv = (Material.alpha - Material.n)/Material.Ks + Material.n/Material.Kf;

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
        ne = 100;
        % column size [m]
        L = 6;
        %%%% solid displacement field
        typeU = 'L3';
        MeshU = Build1DMesh(nsd, ne, L, typeU);
        %%%% fluid pressure field
        typeP = 'L2';
        MeshP = Build1DMesh(nsd, ne, L, typeP);
        %%%% porosity field
        if ~Control.Biotmodel
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
        meshFileNameU = 'Mesh Files\Column2DQ9_Steady_Boone.msh';
        MeshU = BuildMesh_GMSH(meshFileNameU, fieldU, nsd, config_dir, progress_on);
        %%%% pressure field
        fieldP = 'p';
        meshFileNameP = 'Mesh Files\Column2DQ4_Steady_Boone.msh';
        MeshP = BuildMesh_GMSH(meshFileNameP, fieldP, nsd, config_dir, progress_on);
        %%%% porosity field
        if ~Control.Biotmodel
            fieldN = 'n';
            meshFileNameN = 'Mesh Files\Column2DQ4_Steady_Boone.msh';
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
% ---------- test 1: all boundary displacements fixed (represent rigid body
% motion)
% BC.fixed_u_dof1 = MeshU.left_dof;
% BC.fixed_u_dof2 = MeshU.right_dof;
% BC.fixed_u_dof3 = MeshU.bottom_dof;
% BC.fixed_u_dof4 = MeshU.top_dof;
% BC.fixed_u = unique([BC.fixed_u_dof1;BC.fixed_u_dof2;BC.fixed_u_dof3;BC.fixed_u_dof4]);
% % fixed DOF values
% BC.fixed_u_value = zeros(length(BC.fixed_u),1);
% BC.fixed_u_value(1:2:end) = 0;
% BC.fixed_u_value(2:2:end) = 1;

% ---------- test 2: linear displacement in x (represent constant strain
% state)
% BC.ux = @(x) x;
% BC.fixed_u1 = MeshU.xdofs.';
% BC.fixed_u_value1 = BC.ux(MeshU.coords(:,1));
% BC.fixed_u2 = MeshU.ydofs.';
% BC.fixed_u_value2 = zeros(length(BC.fixed_u2),1);
% BC.fixed_u = [BC.fixed_u1; BC.fixed_u2];
% BC.fixed_u_value = [BC.fixed_u_value1; BC.fixed_u_value2];

% ---------- test 3: linear displacement in y (represent constant strain
% state)
% BC.uy = @(y) y;
% BC.fixed_u1 = MeshU.xdofs.';
% BC.fixed_u_value1 = zeros(length(BC.fixed_u1),1);
% BC.fixed_u2 = MeshU.ydofs.';
% BC.fixed_u_value2 = BC.uy(MeshU.coords(:,2));
% BC.fixed_u = [BC.fixed_u1; BC.fixed_u2];
% BC.fixed_u_value = [BC.fixed_u_value1; BC.fixed_u_value2];

% ---------- test 4: fix bottom y, left corner x and y
% bottomleftnode  = find(MeshU.coords(MeshU.bottom_nodes,1) == min(MeshU.coords(:,1)));
% BC.fixed_u = [MeshU.bottom_nodes*2-1; bottomleftnode*2];
% BC.fixed_u_value = zeros(length(BC.fixed_u),1);

% ---------- test 5: fix left x, left corner x and y
% bottomleftnode  = find(MeshU.coords(MeshU.bottom_nodes,1) == min(MeshU.coords(:,1)));
% BC.fixed_u = [MeshU.left_nodes*2-1; bottomleftnode*2];
% BC.fixed_u_value = zeros(length(BC.fixed_u),1);

% ---------- test 6: fix left corner x and y, left x and bottom y
% bottomleftnode  = find(MeshU.coords(MeshU.bottom_nodes,1) == min(MeshU.coords(:,1)));
% BC.fixed_u = [MeshU.left_nodes*2-1; MeshU.bottom_nodes*2; bottomleftnode*2];
% BC.fixed_u_value = zeros(length(BC.fixed_u),1);

% ---------- test 7: fix left and right x, bottom y
BC.fixed_u = [MeshU.left_nodes*2-1; MeshU.right_nodes*2-1; MeshU.bottom_nodes*2];
BC.fixed_u_value = zeros(length(BC.fixed_u),1);

% ---------- test 8: fix bottom
% BC.fixed_u = [MeshU.bottom_nodes*2-1; MeshU.bottom_nodes*2];
% BC.fixed_u_value = zeros(length(BC.fixed_u),1);


% ---------- end of testing -----------------------------------------------

% free nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid
BC.p = @(x) x;
% pressure fixed at all nodes
BC.fixed_p = (1:MeshP.nDOF);
% fixed DOF values
% BC.fixed_p_value = zeros(length(BC.fixed_p),1);
BC.fixed_p_value = BC.p(MeshP.coords(:,2));
% free nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% BC.tractionNodes = [];
% traction interpolation (needed for traction applied in wells); 1 - true, 0 - false
BC.tractionInterp = 0;

% prescribed traction [GN/m2]
BC.traction = 1e-5;

% ---------- test 1: load at top nodes
BC.tractionNodes = MeshU.top_nodes;
Force = BC.traction * max(MeshU.coords(:,1))/((length(MeshU.top_nodes) - 1)/2);
BC.tractionForce = zeros(length(BC.tractionNodes),2);
% Q9 elements for displacement field
for n = 1:length(BC.tractionForce)
    if any(BC.tractionNodes(n) == MeshU.conn(:,1:4),'all') % then node is a corner node
        BC.tractionForce(n,:) = [0, Force/3];
    else % then node is a midside node
        BC.tractionForce(n,:) = [0, Force*2/3];
    end
end
% find the nodes in the left and right top corners
lefttopnode = find(MeshU.coords(BC.tractionNodes,1) == min(MeshU.coords(:,1)));
righttopnode  = find(MeshU.coords(BC.tractionNodes,1) == max(MeshU.coords(:,1)));
% correcting force at nodes
BC.tractionForce(lefttopnode,2) = BC.tractionForce(lefttopnode,2)/2;
BC.tractionForce(righttopnode,2) = BC.tractionForce(righttopnode,2)/2;

% ---------- test 2: load at right nodes
% BC.tractionNodes = MeshU.right_nodes;
% Force = BC.traction * max(MeshU.coords(:,2))/((length(MeshU.right_nodes) - 1)/2);
% BC.tractionForce = zeros(length(BC.tractionNodes),2);
% % Q9 elements for displacement field
% for n = 1:length(BC.tractionForce)
%     if any(BC.tractionNodes(n) == MeshU.conn(:,1:4),'all') % then node is a corner node
%         BC.tractionForce(n,:) = [Force/3, 0];
%     else % then node is a midside node
%         BC.tractionForce(n,:) = [Force*2/3, 0];
%     end
% end
% % find the nodes in the top and bottom right corners
% rightbottomnode = find(MeshU.coords(BC.tractionNodes,2) == min(MeshU.coords(:,2)));
% righttopnode  = find(MeshU.coords(BC.tractionNodes,2) == max(MeshU.coords(:,2)));
% % correcting force at nodes
% BC.tractionForce(rightbottomnode,1) = BC.tractionForce(rightbottomnode,1)/2;
% BC.tractionForce(righttopnode,1) = BC.tractionForce(righttopnode,1)/2;

% ---------- end of testing -----------------------------------------------

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

%% Problem type
% 1 = quasi-steady/transient problem (no acceleration and pressure change)
% 0 = dynamic problem (acceleration/intertia terms included)
Control.steady = 1;

%% Solution parameters
Control.dt = 1e-2;  % time step
Control.tend = 10;   % final simulation time

Control.beta = 1; % beta-method time discretization -- beta = 1 Backward Euler; beta = 0.5 Crank-Nicolson

Control.plotu = 1; % dof y of node 63 (x = 0.05m, y = 3m)
Control.plotp = 1; % dof of node 35 (x = 0.05m, y = 3m)

% plot analytical solution (valid for 1D problems with Material.Minv == 0)
Control.plotansol = 0; % 1 = true; 0 = false

% solve in the frequency domain
Control.freqDomain = 0;  % 1 = true; 0 = false

end