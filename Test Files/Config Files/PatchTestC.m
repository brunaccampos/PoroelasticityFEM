function [Material, MeshU, MeshP, BC, Control] = PatchTestC(config_dir, progress_on, meshfilename)
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

%% Poroelasticity
% porous media permeability [m2/Pa s]
Material.kf = 0;
% 1/Q (related to storage coefficient)
Material.Qinv = 0;
% Biot's coefficient
Material.alpha = 0;

%% Material properties
% elasticity modulus [Pa]
Material.E = 2540;
% Poisson's ratio
Material.nu = 0.3;

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
        
        % solid displacement field
        typeU = 'L3';
        MeshU = Build1DMesh(nsd, ne, L, typeU);
        
        % fluid pressure field
        typeP = 'L2';
        MeshP = Build1DMesh(nsd, ne, L, typeP);
    case 'Gmsh'
        % Version 2 ASCII
        % number of space dimensions
        nsd = 2;
        
        % field
        fieldU = 'u';
        % build mesh displacement field
        meshFileNameU = meshfilename;
        MeshU = BuildMesh_GMSH(meshFileNameU, fieldU, nsd, config_dir, progress_on);
        
        % field
        fieldP = 'p';
        % build mesh pressure field
        meshFileNameP = 'Mesh Files\PatchTest.msh';
        MeshP = BuildMesh_GMSH(meshFileNameP, fieldP, nsd, config_dir, progress_on);
end

%% Dirichlet BCs
% traction [Pa] (same in both directions)
BC.traction = 3.495;
% displacements according to exact solution
BC.ux = @(x) (1-Material.nu)*BC.traction/Material.E*x(:,1);
BC.uy = @(x) (1-Material.nu)*BC.traction/Material.E*x(:,2);

BC.fixed_u_dof1 = [MeshU.left_nodes*2-1];
BC.fixed_u_dof2 = MeshU.bottom_nodes*2;
BC.fixed_u = [BC.fixed_u_dof1; BC.fixed_u_dof2];

% prescribed displacement
BC.fixed_u_value = zeros(length(BC.fixed_u),1);
BC.fixed_u_value1 = BC.ux([MeshU.coords(MeshU.left_nodes,1),MeshU.coords(MeshU.left_nodes,2)]);
BC.fixed_u_value2 = BC.uy([MeshU.coords(MeshU.bottom_nodes,1),MeshU.coords(MeshU.bottom_nodes,2)]);
BC.fixed_u_value = [BC.fixed_u_value1; BC.fixed_u_value2];

% prescribed pressure
BC.fixed_p = 1:MeshP.nDOF;
BC.fixed_p_value = zeros(length(BC.fixed_p),1);
BC.fixed_p_value = zeros(length(BC.fixed_p),1);

%% Neumann BCs
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

% prescribed flux
BC.fixed_fp = 1:MeshP.nDOF;

% free nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);
BC.free_fp = setdiff(MeshP.DOF, BC.fixed_fp);

% point load [N]
BC.pointLoad = [];

% body force
BC.b = @(x)[];

%% Quadrature order
Control.nq = 2;

%% Problem type
% 1 = quasi-steady state problem (no solid velocity, acceleration, and pressure
% change)
% 0 = transient problem (velocity and acceleration included)
Control.steady = 1;

%% Solution parameters
Control.dt = 1;  % time step

end