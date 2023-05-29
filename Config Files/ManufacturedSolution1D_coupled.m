function [Material, MeshU, MeshP, MeshN, BC, Control] = ManufacturedSolution1D_coupled(config_dir, progress_on)
% ------------------------------------------------------------------------
% Manufactured solution for element mesh size convergence study
% u = sin(x)
% p = sin(x)
% ------------------------------------------------------------------------
% Adapted from https://github.com/GCMLab (Acknowledgements: Bruce Gee)
% ------------------------------------------------------------------------

%% Poroelasticity
% 1/Q (related to storage coefficient)
Material.Minv = 0;
% Biot's coefficient
Material.alpha = 0;
% poroelasticity model
Control.PMmodel = 'Tr1_Biot_UP';

% thickness 
% 1D: cross sectional area [m2]
% 2D: out of plane thickness [m]
Material.t = 1;

% constititive law - 'PlaneStress' or 'PlaneStrain'
% Note: use 'PlaneStrain' for 1D or 2D poroelasticity
Material.constLaw = 'PlaneStress';

%% Material properties
% porous media permeability [m2/Pa s]
Material.kf = 1;
% elasticity modulus [Pa]
Material.E = 1;
% Poisson's ratio
Material.nu = 0.3;

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
        % number of space dimensions
        nsd = 1;
        % number of elements
        ne = 128;
        % column size [m]
        L = 1;
        
        % solid displacement field
        typeU = 'L3';
        fieldU = 'u';
        MeshU = Build1DMesh(nsd, ne, L, typeU, fieldU);
        
        % fluid pressure field
        typeP = 'L2';
        fieldP = 'p';
        MeshP = Build1DMesh(nsd, ne, L, typeP, fieldP);
        
        MeshN = [];
        
    case 'Gmsh'
        % Version 2 ASCII
        % number of space dimensions
        nsd = 2;
        %%%% displacement field
        fieldU = 'u';
        % build mesh displacement field
        meshFileNameU = 'Mesh Files\Manufactured_finerQ4.msh';
        MeshU = BuildMesh_GMSH(meshFileNameU, fieldU, nsd, config_dir, progress_on);
        %%%% pressure field
        fieldP = 'p';
        % build mesh pressure field
        meshFileNameP = 'Mesh Files\Manufactured_finerQ4.msh';
        MeshP = BuildMesh_GMSH(meshFileNameP, fieldP, nsd, config_dir, progress_on);
        %%%% porosity field
        MeshN = [];
end

%% Initial conditions
% displacement
BC.initU = [];

% pressure
BC.initP = [];

%% Find nodes for prescribed BCs
% find top and bottom nodes for displacement field
BC.top_node_u = find(MeshU.coords == max(MeshU.coords));
BC.bottom_node_u = find(MeshU.coords == min(MeshU.coords));

% find top and bottom nodes for pressure field
BC.top_node_p = find(MeshP.coords == max(MeshP.coords));
BC.bottom_node_p = find(MeshP.coords == min(MeshP.coords));

%% Dirichlet BCs - solid
BC.ux = @(x) sin(x);
% column vector of prescribed displacement dof
BC.fixed_u_dof1 = BC.top_node_u;
BC.fixed_u_dof2 = BC.bottom_node_u;
BC.fixed_u = [BC.fixed_u_dof1; BC.fixed_u_dof2];
% prescribed displacement
BC.fixed_u_value = zeros(length(BC.fixed_u),1);
BC.fixed_u_value = BC.ux(MeshU.coords(BC.fixed_u));
% free displacement nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid
BC.p = @(x) sin(x);
% column vector of prescribed displacement dof
BC.fixed_p_dof1 = BC.top_node_p;
BC.fixed_p_dof2 = BC.bottom_node_p;
BC.fixed_p = [BC.fixed_p_dof1; BC.fixed_p_dof2];
% prescribed pressure
BC.fixed_p_value = zeros(length(BC.fixed_p),1);
BC.fixed_p_value = BC.p(MeshP.coords(BC.fixed_p));
% free pressure nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% column vector of prescribed traction nodes
BC.tractionNodes = [];

% body force
BC.b = @(x) sin(x);

% point load [N]
BC.pointLoad = [];

%% Neumann BCs - fluid
% point flux [m/s]
BC.pointFlux = [];

% distributed flux [m/s]
BC.fluxNodes = [];

% flux source
BC.s = @(x) sin(x);

%% Quadrature order
Control.nqU = 3;
Control.nqP = 3;

%% Solution parameters
% tag used for computing analytical solution
% 1 = uncoupled problem (elasticity, heat transfer, etc)
% 0 = coupled problem (Biot, Spanos model)
Control.uncoupled = 1; 

% analytical solution
Control.pan_symb = @(x) sin(x);
Control.p_an = Control.pan_symb(MeshP.coords);
Control.uan_symb = @(x) sin(x);
Control.u_an = Control.uan_symb(MeshU.coords);

Control.dt = 1;  % time step
Control.tend = 1;
Control.step = 1;

Control.beta = 1; % beta-method time discretization -- beta = 1 Backward Euler; beta = 0.5 Crank-Nicolson

Control.plotu = round(MeshU.nn/2);
Control.plotp = round(MeshP.nn/2);

% plot analytical solution (valid for 1D problems with Material.Minv == 0)
Control.plotansol = 1; % 1 = true; 0 = false

% solve in the frequency domain
Control.freqDomain = 0;  % 1 = true; 0 = false

end
