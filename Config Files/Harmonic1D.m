function [Material, MeshU, MeshP, MeshN, BC, Control] = Harmonic1D(config_dir, progress_on)
% Wave propagation in 1D
% Configuration File
% Based on Zienkiewicz (1982) model
% ------------------------------------------------------------------------
% Assumptions/conventions:
% - stress is positive for tension
% - boundary condition for force is based on total stress
% - only solid acceleration is considered (undrained condition; no motions
% of the fluid relative to the solid skeleton can occur)
% - solid grains and fluid are incompressible
% ------------------------------------------------------------------------

%% Poroelasticity model
% 1 - Biot theory
% 0 - Spanos theory (additional porosity equation)
Control.Biotmodel = 1;

%% Material properties - Komijani (2019)
% elasticity modulus [GPa]
Material.E = 14.516e-3;
% Poisson's ratio
Material.nu = 0.3;
% porous media permeability [m2/GPa s]
Material.kf = 1.0194e3;
% dynamic viscosity [GPa s]
Material.mu = 1e-12;
% intrinsic permeability [m2]
Material.k = Material.kf * Material.mu;
% fluid bulk modulus [GPa]
Material.Kf = 2.1;
% solid bulk modulus [GPa]
Material.Ks = 1e11;
% material porosity
Material.n = 0.3;
% Biot's coefficient
Material.alpha = 1;
% fluid density [10^9 kg/m3]
Material.rho_f = 1000e-9;
% solid density [10^9 kg/m3]
Material.rho_s = 2000e-9;
% average density of the medium
Material.rho = Material.n*Material.rho_f + (1-Material.n)*Material.rho_s;
% 1/Q (related to storage coefficient)
Material.Minv = (Material.alpha - Material.n)/Material.Ks + Material.n/Material.Kf;

% lumped mass matrix - 0: false, 1: true
Material.lumpedMass = 0;

% constititive law - 'PlaneStress' or 'PlaneStrain'
% Note: use 'PlaneStrain' for 1D or 2D poroelasticity
Material.constLaw = 'PlaneStrain';

%% Material parameters - harmonic oscillation formulation
% circular frequency
Control.omega = 1; % [Hz]

% Biot (1956)
Material.rho12 = 0; % added mass
Material.rho11 = (1 - Material.n) * Material.rho_s - Material.rho12;
Material.rho22 = Material.n * Material.rho_f - Material.rho12;
Material.b = Material.mu * Material.n^2 / Material.k;

% Atalla (1998)
Material.rho11_tilde = Material.rho11 + Material.b/(1i * Control.omega);
Material.rho22_tilde = Material.rho22 + Material.b/(1i * Control.omega);
Material.rho12_tilde = Material.rho12 - Material.b/(1i * Control.omega);

Material.rho_tilde = Material.rho11_tilde - Material.rho12_tilde^2 / Material.rho22_tilde;
Material.M = 1/Material.Minv;
Material.R = Material.n * Material.Kf * (1 - Material.M * (Material.alpha - Material.n) / Material.Ks);
Material.gamma_tilde = Material.n * (Material.rho12_tilde / Material.rho22_tilde - Material.M/Material.R);

%% In situ stress field
% [GPa]
BC.S0 = [];

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
        ne = 100;
        % column size [m]
        L = 10;
        %%%% solid displacement field
        typeU = 'L3';
        fieldU = 'u';
        MeshU = Build1DMesh(nsd, ne, L, typeU, fieldU);
        %%%% fluid pressure field
        typeP = 'L2';
        fieldP = 'p';
        MeshP = Build1DMesh(nsd, ne, L, typeP, fieldP);
        %%%% porosity field
        if ~Control.Biotmodel
            typeN = 'L2';
            fieldN = 'n';
            MeshN = Build1DMesh(nsd, ne, L, typeN, fieldN);
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
% displacement u=0 at the bottom
BC.fixed_u = (BC.bottom_node_u);
BC.fixed_u_value = zeros(length(BC.fixed_u),1);
% free displacement nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid
%   pressure p=0 at the top
BC.fixed_p = (BC.top_node_p);
BC.fixed_p_value = zeros(length(BC.fixed_p),1);
% free pressure nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% point load [GN]
BC.pointLoadValue = -3000e-9;
BC.pointLoadNodes = BC.top_node_u;
BC.pointLoad = zeros(MeshU.nDOF,1);
BC.pointLoad(BC.pointLoadNodes) = BC.pointLoadValue;

% distributed load [GN/m2]
BC.tractionNodes = [];

% body force [GN/m3]
BC.b = @(x)[];  

%% Neumann BCs - fluid
% point flux [m/s]
BC.pointFluxValue = 0;
BC.pointFluxNodes = BC.bottom_node_p;
BC.pointFlux = zeros(MeshP.nDOF,1);
BC.pointFlux(BC.pointFluxNodes) = BC.pointFluxValue;

% distributed flux [m3/s]
BC.fluxNodes = [];

% flux source [m3/s/m3]
BC.s = @(x)[]; 

%% Quadrature order
Control.nqU = 2;
Control.nqP = 2;

%% Problem type
% 1 = quasi-steady/transient problem (no acceleration and pressure change)
% 0 = dynamic problem (acceleration/intertia terms included)
Control.steady = 0;

% tag used for computing analytical solution
% 1 = uncoupled problem (elasticity, heat transfer, etc)
% 0 = coupled problem (Biot, Spanos model)
Control.uncoupled = 0; 

%% Solution parameters
Control.dt = 1e-2;  % time step
Control.tend = 10;   % final simulation time

Control.plotu = find(MeshU.coords == 5); % x = 5m
Control.plotp = find(MeshP.coords == 5); % x = 5m

% plot analytical solution (valid for 1D problems with Material.Minv == 0)
Control.plotansol = 0; % 1 = true; 0 = false

% solve in the frequency domain
Control.freqDomain = 0;  % 1 = true; 0 = false

%% Time discretization parameters
% Newmark method
Control.beta = 0.7;
Control.gamma = 0.7;
Control.theta = 0.7;

end