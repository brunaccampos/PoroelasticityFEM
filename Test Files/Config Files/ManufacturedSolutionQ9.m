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

%% Material properties
% elasticity modulus [Pa]
Material.E = 2230;
% Poisson's ratio
Material.nu = 0.3;
% porous media permeability [m2/Pa s]
Material.kf = 0;
% 1/Q (related to storage coefficient)
Material.Minv = 0;
% Biot's coefficient
Material.alpha = 0;

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
% auxiliar vector
displ = zeros(length(BC.fixed_u),1);
displ(1:2:end) = BC.ux([MeshU.coords(BC.fixed_u(2:2:end)/2,1),MeshU.coords(BC.fixed_u(2:2:end)/2,2)]);
displ(2:2:end) = BC.uy([MeshU.coords(BC.fixed_u(2:2:end)/2,1),MeshU.coords(BC.fixed_u(2:2:end)/2,2)]);
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
% column vector of prescribed traction nodes
BC.tractionNodes = MeshU.right_nodes;
% prescribed traction [t1x t1y;t2x t2y;...] [N]
Fnode = 1/(length(BC.tractionNodes) - 1);
BC.tractionForce = Fnode*[zeros(size(BC.tractionNodes)), zeros(size(BC.tractionNodes))];

% time dependent vector
BC.tractionForce = @(t) BC.tractionForce;

% body force
E = Material.E;
nu = Material.nu;
BC.b = @(x,t)[-E / (1-nu^2)  * ( 20*x(1).^3 + 3*nu*x(2).^2              + (1-nu)/2*( 6*x(1).*x(2) - 30*x(2).^4 + 3*x(2).^2));
    -E / (1-nu^2)  * ( (1-nu)/2*( 3*x(2).^2  + 20*x(1).^3)    + 3*nu*x(2).^2 + 6*x(1).*x(2) - 30*x(2).^4 )];

% point load [N]
BC.pointLoad = [];

%% Neumann BCs - fluid
% point flux [m/s]
BC.pointFlux = [];

% distributed flux [m/s]
BC.fluxNodes = [];

% flux source
BC.s = @(x,t)[]; 

%% Quadrature order
Control.nqU = 4;
Control.nqP = 4;

%% Time step controls
Control.dt = 1;  % time step
Control.tend = 1; % final simulation time

% Beta method
% beta = 1 Backward Euler; beta = 0.5 Crank-Nicolson
Control.beta = 1; 

end
