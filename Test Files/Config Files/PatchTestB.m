function [Material, MeshU, MeshP, MeshN, BC, Control] = PatchTestB(config_dir, progress_on, meshfilename)
% ------------------------------------------------------------------------
% For Patch Test B, only nodes 1-8 (nodes in the boundaries) are restrained
% with their displacements specified according to the exact solution. The error between the FEA
% and exact solutions is then calculated. The FEA approximate solution
% should be exact.
% ------------------------------------------------------------------------
% Adapted from https://github.com/GCMLab (Acknowledgements: Bruce Gee)
% ------------------------------------------------------------------------

%% Poroelasticity
% porous media permeability [m2/Pa s]
Material.kf = 0;
% 1/Q (related to storage coefficient)
Material.Minv = 0;
% Biot's coefficient
Material.alpha = 0;
% poroelasticity model
Control.PMmodel = 'Tr1_Biot_UP';

%% Material properties
% elasticity modulus [Pa]
Material.E = 2540;
% Poisson's ratio
Material.nu = 0.3;

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

%% Dirichlet BCs - solid
% displacements according to exact solution
BC.ux = @(x) (1-Material.nu)*BC.traction/Material.E*x(:,1);
BC.uy = @(x) (1-Material.nu)*BC.traction/Material.E*x(:,2);
% column vector of prescribed displacement dof
BC.fixed_u_dof1 = MeshU.left_dof;
BC.fixed_u_dof2 = MeshU.right_dof;
BC.fixed_u_dof3 = MeshU.bottom_dof;
BC.fixed_u_dof4 = MeshU.top_dof;
BC.fixed_u = unique([BC.fixed_u_dof1;BC.fixed_u_dof2;BC.fixed_u_dof3;BC.fixed_u_dof4]);
% prescribed displacement
BC.fixed_u_value = zeros(length(BC.fixed_u),1);
% auxiliar vector
displ = zeros(length(BC.fixed_u),1);
displ(1:2:end) = BC.ux([MeshU.coords(BC.fixed_u(2:2:end)/2,1), MeshU.coords(BC.fixed_u(2:2:end)/2,2)]);
displ(2:2:end) = BC.uy([MeshU.coords(BC.fixed_u(2:2:end)/2,1), MeshU.coords(BC.fixed_u(2:2:end)/2,2)]);
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
% traction interpolation (needed for traction applied in wells); 1 - true, 0 - false
BC.tractionInterp = 0;

BC.tractionNodes = MeshU.right_nodes;
Force = 1/(length(BC.tractionNodes) - 1);
BC.tractionForce = Force*[zeros(size(BC.tractionNodes)), zeros(size(BC.tractionNodes))];
% find the nodes in the top right and bottom right corners
toprightnode = find(MeshU.coords(BC.tractionNodes,2) == max(MeshU.coords(:,2)));
botrightnode = find(MeshU.coords(BC.tractionNodes,2) == min(MeshU.coords(:,2)));

BC.tractionForce(toprightnode,1) = BC.tractionForce(toprightnode,1)/2;
BC.tractionForce(botrightnode,1) = BC.tractionForce(botrightnode,1)/2;

% point load [N]
BC.pointLoad = [];

% body force
BC.b = @(x,t)[];

%% Neumann BCs - fluid
% point flux [m/s]
BC.pointFlux = [];

% distributed flux [m/s]
BC.fluxNodes = [];

% flux source
BC.s = @(x,t)[]; 

%% Quadrature order
Control.nqU = 2;
Control.nqP = 2;

%% Problem type
% 1 = quasi-steady/transient problem (no acceleration and pressure change)
% 0 = dynamic problem (acceleration/intertia terms included)
Control.steady = 1;

%% Solution parameters
Control.dt = 1;  % time step
Control.t = 0; % time variable
Control.step = 1; % total simulation time

Control.beta = 1; % beta-method time discretization -- beta = 1 Backward Euler; beta = 0.5 Crank-Nicolson

end