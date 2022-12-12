function [Material, MeshU, MeshP, BC, Control] = PatchTestA(config_dir, progress_on, meshfilename)
% ------------------------------------------------------------------------
% Patch Test A: all nodes are restrained and nodal displacement values
% are specfied according to the exact solution. The error between the FEA
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
        meshFileNameP = meshfilename;
        MeshP = BuildMesh_GMSH(meshFileNameP, fieldP, nsd, config_dir, progress_on);
end

%% Dirichlet BCs
% traction [Pa] (same in both directions)
BC.traction = 3.495;
% displacements according to exact solution
BC.ux = @(x) (1-Material.nu)*BC.traction/Material.E*x(:,1);
BC.uy = @(x) (1-Material.nu)*BC.traction/Material.E*x(:,2);
% fixed nodes
BC.fixed_u = 1:MeshU.nDOF;
% prescribed displacements
BC.fixed_u_value = zeros(length(BC.fixed_u),1);
BC.fixed_u_value(1:2:end) = BC.ux(MeshU.coords);
BC.fixed_u_value(2:2:end) = BC.uy(MeshU.coords);
% prescribed pressure
BC.fixed_p = 1:MeshP.nDOF;
BC.fixed_p_value = zeros(length(BC.fixed_p),1);
BC.fixed_p_value = zeros(length(BC.fixed_p),1);

%% Neumann BCs
BC.tractionNodes = MeshU.right_nodes;
Force = 1/(length(BC.tractionNodes) - 1);
BC.tractionForce = Force*[zeros(size(BC.tractionNodes)), zeros(size(BC.tractionNodes))];
% find the nodes in the top right and bottom right corners
toprightnode = find(MeshU.coords(BC.tractionNodes,2) == max(MeshU.coords(:,2)));
botrightnode = find(MeshU.coords(BC.tractionNodes,2) == min(MeshU.coords(:,2)));

BC.tractionForce(toprightnode,1) = BC.tractionForce(toprightnode,1)/2;
BC.tractionForce(botrightnode,1) = BC.tractionForce(botrightnode,1)/2;

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