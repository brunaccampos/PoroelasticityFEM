function [Material, MeshU, MeshP, MeshN, BC, Control] = PatchTestC(config_dir, progress_on, ~, ~)
% ------------------------------------------------------------------------
% Patch Test C is performed with node 1 fully restrained and nodes 4 and 8
% restrained only in the x -direction. Nodal forces are applied to nodes 2,
% 3, and 6 in accordance with the values generated through the boundary
% tractions by sigma(x)=2. The error between the FEA
% and exact solutions is then calculated. The FEA approximate solution
% should be exact.
% ------------------------------------------------------------------------
% Adapted from https://github.com/GCMLab (Acknowledgements: Bruce Gee)
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
Control.PMmodel = 'Tr_BT_UP';

%% Material properties
% elasticity modulus [Pa]
Material.M(1).E = 2540;
% Poisson's ratio
Material.M(1).nu = 0.3;
% porous media permeability [m2/Pa s]
Material.M(1).kf = 0;
% 1/Q (related to storage coefficient)
Material.M(1).Minv = 0;
% Biot's coefficient
Material.M(1).alpha = 0;

% thickness 
% 1D: cross sectional area [m2]
% 2D: out of plane thickness [m]
Material.t = 1;

% constititive law - 'PlaneStress' or 'PlaneStrain'
Material.constLaw = 'PlaneStress';

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
        L = 6;
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
        % build mesh displacement field
        meshFileNameU = 'Mesh Files\PatchTest.msh';
        MeshU = BuildMesh_GMSH(meshFileNameU, fieldU, nsd, config_dir, progress_on);
        %%%% pressure field
        fieldP = 'p';
        % build mesh pressure field
        meshFileNameP = 'Mesh Files\PatchTest.msh';
        MeshP = BuildMesh_GMSH(meshFileNameP, fieldP, nsd, config_dir, progress_on);
        %%%% porosity field
        if contains(Control.PMmodel, 'UPN')
            fieldN = 'n';
            meshFileNameN = 'Mesh Files\PatchTest.msh';
            MeshN = BuildMesh_GMSH(meshFileNameN, fieldN, nsd, config_dir, progress_on);
        else
            MeshN = [];
        end
end

% traction [Pa] (same in both directions)
BC.traction = 3.495;

%% Dirichlet BCs
% displacements according to exact solution
BC.ux = @(x) (1-Material.M(1).nu)*BC.traction/Material.M(1).E*x(:,1);
BC.uy = @(x) (1-Material.M(1).nu)*BC.traction/Material.M(1).E*x(:,2);
% fixed nodes
BC.fixed_u_dof1 = MeshU.left_nodes*2-1;
BC.fixed_u_dof2 = MeshU.bottom_nodes*2;
BC.fixed_u = [BC.fixed_u_dof1; BC.fixed_u_dof2];
% prescribed displacement
BC.fixed_u_value = zeros(length(BC.fixed_u),1);
% auxiliar vector
displ1 = BC.ux([MeshU.coords(MeshU.left_nodes,1),MeshU.coords(MeshU.left_nodes,2)]);
displ2 = BC.uy([MeshU.coords(MeshU.bottom_nodes,1),MeshU.coords(MeshU.bottom_nodes,2)]);
displ = [displ1; displ2];
% atribute vector to time dependent BC
BC.fixed_u_value = @(t) displ;
% free displacement nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid
% prescribed pressure
BC.fixed_p = 1:MeshP.nDOF;
BC.fixed_p_value = @(t) zeros(length(BC.fixed_p),1);
% free pressure nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% column vector of prescribed traction nodes
toprightnode = MeshU.right_nodes(MeshU.coords(MeshU.right_nodes,2) == max(MeshU.coords(:,2)));
index_right = MeshU.right_nodes ~= toprightnode;
index_top   = MeshU.top_nodes   ~= toprightnode;
BC.tractionNodes = [MeshU.right_nodes(index_right);  MeshU.top_nodes(index_top); toprightnode];

% prescribed traction
Fright = BC.traction * max(MeshU.coords(:,2))/(length(MeshU.right_nodes) - 1);
Ftop   = BC.traction * max(MeshU.coords(:,1))/(length(MeshU.top_nodes)   - 1);
BC.tractionForce = [Fright*ones(size(MeshU.right_nodes(index_right))), zeros(size(MeshU.right_nodes(index_right))); % right side nodes
    zeros(size(MeshU.top_nodes(index_top))), Ftop*ones(size(MeshU.top_nodes(index_top))); % top side nodes
    Fright*1/2, Ftop*1/2]; % top right node

% find the nodes in the top left and bottom right corners
botrightnode = find(MeshU.coords(BC.tractionNodes,2) == min(MeshU.coords(:,2)));
topleftnode  = find(MeshU.coords(BC.tractionNodes,1) == min(MeshU.coords(:,1)));

BC.tractionForce(botrightnode,1) = BC.tractionForce(botrightnode,1)/2;
BC.tractionForce(topleftnode,2) = BC.tractionForce(topleftnode,2)/2;

% point load [N]
BC.pointLoad = @(t)[];

% body force
BC.b = @(x,t)[];

%% Neumann BCs - fluid
% point flux [m/s]
BC.pointFlux = @(t)[];

% distributed flux [m/s]
BC.fluxNodes = [];

% flux source
BC.s = @(x,t)[]; 

%% Quadrature order
Control.nqU = 2;
Control.nqP = 2;

%% Time step controls
Control.dt = 1;  % time step
Control.tend = 1; % final simulation time

% Beta method
% beta = 1 Backward Euler; beta = 0.5 Crank-Nicolson
Control.beta = 1; 

%% Plot data
% DOF to plot graphs
Control.plotu = 1;
Control.plotp = 1;

end