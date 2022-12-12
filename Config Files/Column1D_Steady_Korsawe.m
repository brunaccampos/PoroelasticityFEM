function [Material, MeshU, MeshP, MeshN, BC, Control] = Column1D_Steady_Korsawe(config_dir, progress_on)
% Column Consolidation 1D simulation
% Configuration File
% Based on Korsawe (2006) model
% ------------------------------------------------------------------------
% Assumptions/conventions:
% - stress is positive for tension
% - boundary condition for force is based on total stress
% - no acceleration terms for solid or fluid
% - solid velocity is neglected
% - fluid and solid grains are incompressible
% - porosity is constant in space and varies over time
% ------------------------------------------------------------------------

%% Poroelasticity model
% 1 - Biot theory
% 0 - Spanos theory (additional porosity equation)
Control.Biotmodel = 0;

%% Material properties - Korsawe (2006)
% elasticity modulus [Pa]
Material.E = 3e4;
% Poisson's ratio
Material.nu = 0.2;
% intrinsic permeability [m2]
Material.k = 1e-10;
% dynamic viscosity [Pa s]
Material.mu = 1e-3;
% porous media permeability [m2/Pa s]
Material.kf = Material.k/Material.mu;
% Biot's coefficient
Material.alpha = 1;
% 1/Q (related to storage coefficient)
Material.Qinv = 0;

%%% Aditional parameters from Komijani (2019) for varying porosity model
% fluid bulk modulus [Pa]
Material.Kf = 2.1e9;
% solid bulk modulus [Pa]
Material.Ks = 1e20;
% material porosity
Material.n = 0.3;

%% Spanos material parameters
n = Material.n;
% % porosity coefficients from Spanos (1989)
Material.deltaF = n*(1-n)/(Material.Ks*(n/Material.Kf + (1-n)/Material.Ks));
Material.deltaS = n*(1-n)/(Material.Kf*(n/Material.Kf + (1-n)/Material.Ks));
% Material.deltaF = 0;
% Material.deltaS = 0;

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
        ne = 10;
        % column size [m]
        L = 1;
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
        meshFileNameU = 'Column2DQ9.msh';
        MeshU = BuildMesh_GMSH(meshFileNameU, fieldU, nsd, config_dir, progress_on);
        %%%% pressure field
        fieldP = 'p';
        meshFileNameP = 'Column2DQ4.msh';
        MeshP = BuildMesh_GMSH(meshFileNameP, fieldP, nsd, config_dir, progress_on);
        %%%% porosity field
        if ~Control.Biotmodel
            fieldN = 'n';
            meshFileNameN = 'Column2DQ4.msh';
            MeshN = BuildMesh_GMSH(meshFileNameN, fieldN, nsd, config_dir, progress_on);
        else
            MeshN = [];
        end
end

%% Boundary conditions

% find top and bottom nodes for displacement field
BC.top_node_u = find(MeshU.coords == min(MeshU.coords));
BC.bottom_node_u = find(MeshU.coords == max(MeshU.coords));

% find top and bottom nodes for pressure field
BC.top_node_p = find(MeshP.coords == min(MeshP.coords));
BC.bottom_node_p = find(MeshP.coords == max(MeshP.coords));

% fixed nodes
%   displacement u=0 at the bottom
%   pressure p=0 at the top
%   no pressure gradient fp=0 at the bottom
BC.fixed_u = (BC.bottom_node_u);
BC.fixed_p = (BC.top_node_p);
BC.fixed_fp = (BC.bottom_node_p);

% free nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);
BC.free_fp = setdiff(MeshP.DOF, BC.fixed_fp);

%% Dirichlet BCs
BC.fixed_u_value = zeros(length(BC.fixed_u),1);
BC.fixed_p_value = zeros(length(BC.fixed_p),1);

%% Neumann BCs
% point load [N]
BC.pointLoadValue = 1e3;
BC.pointLoadNodes = BC.top_node_u;
BC.pointLoad = zeros(MeshU.nDOF,1);
BC.pointLoad(BC.pointLoadNodes) = BC.pointLoadValue;

% distributed load [N/m]
BC.tractionNodes = [];

% body force
BC.b = @(x)[];  

%% Porosity BCs
if ~Control.Biotmodel
    BC.fixed_n = [];
    BC.free_n = setdiff(MeshN.DOF, BC.fixed_n);
    BC.fixed_n_value = zeros(length(BC.fixed_n),1);
end

%% Quadrature order
Control.nq = 2;

%% Problem type
% 1 = quasi-steady state problem (no solid velocity, acceleration, and pressure
% change)
% 0 = transient problem (velocity and acceleration included)
Control.steady = 1;

%% Solution parameters
Control.dt = 1;  % time step
Control.tend = 10;   % final simulation time
Control.plotu = round(length(MeshU.coords)/2);
Control.plotp = round(length(MeshP.coords)/2);

end