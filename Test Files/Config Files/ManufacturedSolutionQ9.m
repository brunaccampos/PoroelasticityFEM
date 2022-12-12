function [Material, MeshU, MeshP, MeshN, BC, Control] = ManufacturedSolutionQ9(config_dir, progress_on, meshfilename)
% ------------------------------------------------------------------------
% Test 4 calculates the convergence rates of a uniform Q4 mesh using a
% manufactured solution in which
% ux = x^5 + x*y^3 - y^6
% uy = x^5 + x*y^3 - y^6
% under plane stress conditions
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
Control.Biotmodel = 1;

%% Material properties
% elasticity modulus [Pa]
Material.E = 2230;
% Poisson's ratio
Material.nu = 0.3;

% constititive law - 'PlaneStress' or 'PlaneStrain'
Material.constLaw = 'PlaneStress';

% lumped mass matrix - 0: false, 1: true
Material.lumpedMass = 0;

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
        
        MeshN = [];
end

%% Dirichlet BCs - solid
BC.ux = @(x) x(:,1).^5 + x(:,1).*x(:,2).^3 - x(:,2).^6;
BC.uy = @(x) x(:,1).^5 + x(:,1).*x(:,2).^3 - x(:,2).^6;
% column vector of prescribed displacement dof
BC.fixed_u_dof1 = MeshU.left_dof;
BC.fixed_u_dof2 = MeshU.right_dof;
BC.fixed_u_dof3 = MeshU.bottom_dof;
BC.fixed_u_dof4 = MeshU.top_dof;
BC.fixed_u = unique([BC.fixed_u_dof1; BC.fixed_u_dof2; BC.fixed_u_dof3; BC.fixed_u_dof4]);
% prescribed displacement
BC.fixed_u_value = zeros(length(BC.fixed_u),1);
BC.fixed_u_value(1:2:end) = BC.ux([MeshU.coords(BC.fixed_u(2:2:end)/2,1),MeshU.coords(BC.fixed_u(2:2:end)/2,2)]);
BC.fixed_u_value(2:2:end) = BC.uy([MeshU.coords(BC.fixed_u(2:2:end)/2,1),MeshU.coords(BC.fixed_u(2:2:end)/2,2)]);
% free displacement nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid
% prescribed pressure
BC.fixed_p = 1:MeshP.nDOF;
BC.fixed_p_value = zeros(length(BC.fixed_p),1);
% free pressure nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% traction interpolation (needed for traction applied in wells); 1 - true, 0 - false
BC.tractionInterp = 0;
% column vector of prescribed traction nodes
BC.tractionNodes = MeshU.right_nodes;
% prescribed traction [t1x t1y;t2x t2y;...] [N]
Fnode = 1/(length(BC.tractionNodes) - 1);
BC.tractionForce = Fnode*[zeros(size(BC.tractionNodes)), zeros(size(BC.tractionNodes))];

% NOTE: point loads at any of the element nodes can also be
% added as a traction.

% magnitude of distributed body force [N/m] [bx;by] according to
% the manufactured solution:
% bx = -E / (1-v^2) * ( 20x^3 + 3vy^2          + (1-v)/2*[ 6xy - 30y^4 + 3y^2 ])
% by = -E/  (1-v^2) * ( (1-v)/2*[3y^2 + 20x^3] +           3vy^2 + 6xy - 30y^4 )
% 1D: [N/m], 2D: [N/m2]
% NOTE: if no body force, use '@(x)[]'
% NOTE: anonymous functions is defined with respect to the
%      variable x,  which is a vector [x(1) x(2)] = [x y]
E = Material.E;
nu = Material.nu;

% body force
BC.b = @(x)[-E / (1-nu^2)  * ( 20*x(1).^3 + 3*nu*x(2).^2              + (1-nu)/2*( 6*x(1).*x(2) - 30*x(2).^4 + 3*x(2).^2));
    -E / (1-nu^2)  * ( (1-nu)/2*( 3*x(2).^2  + 20*x(1).^3)    + 3*nu*x(2).^2 + 6*x(1).*x(2) - 30*x(2).^4 )];

% point load [N]
BC.pointLoad = [];

%% Neumann BCs - fluid
% point flux [m/s]
BC.pointFlux = [];

% distributed flux [m/s]
BC.fluxNodes = [];

% flux source
BC.s = @(x)[]; 

%% Quadrature order
Control.nqU = 4;
Control.nqP = 4;

%% Problem type
% 1 = quasi-steady/transient problem (no acceleration and pressure change)
% 0 = dynamic problem (acceleration/intertia terms included)
Control.steady = 1;

%% Solution parameters
Control.dt = 1;  % time step
Control.step = 1; % total simulation time

Control.beta = 1; % beta-method time discretization -- beta = 1 Backward Euler; beta = 0.5 Crank-Nicolson

end
