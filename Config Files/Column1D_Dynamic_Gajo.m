function [Material, MeshU, MeshP, MeshN, BC, Control] = Column1D_Dynamic_Gajo(config_dir, progress_on)
% Column Consolidation 1D simulation
% Configuration File
% Based on Zienkiewicz (1982) model
%
% Assumptions/conventions:
% - stress is positive for tension
% - boundary condition for force is based on total stress
% - only solid acceleration is considered (undrained condition; no motions
% of the fluid relative to the solid skeleton can occur)
% - solid grains and fluid are incompressible

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
Control.PMmodel = 'Dyn4_Biot_UPU';

%% Material properties - Tasiopoulou (2015)
% elasticity modulus [GPa]
Material.E = 0.2321e2;
% Poisson's ratio
Material.nu = 0.171;
% porous media permeability [m2/GPa s]
Material.kf = 0.1480e4;
% dynamic viscosity [GPa s]
Material.mu = 1e-12;
% intrinsic permeability [m2]
Material.k = Material.kf * Material.mu;
% fluid bulk modulus [GPa]
Material.Kf = 2.2;
% solid bulk modulus [GPa]
Material.Ks = 36;
% material porosity
Material.n = 0.18;
% Biot's coefficient
Material.alpha = 0.677;
% fluid density [10^9 kg/m3]
Material.rho_f = 1000e-9;
% solid density [10^9 kg/m3]
Material.rho_s = 2660e-9;
% average density of the medium
Material.rho = Material.n*Material.rho_f + (1-Material.n)*Material.rho_s;
% 1/Q (related to storage coefficient)
Material.Minv = (Material.alpha - Material.n)/Material.Ks + Material.n/Material.Kf;
% fluid bulk viscosity [GPa s]
Material.xif = 2.8e-12; % (Quiroga-Goode, 2005)

% thickness 
% 1D: cross sectional area [m2]
% 2D: out of plane thickness [m]
Material.t = 1;

% constititive law - 'PlaneStress' or 'PlaneStrain'
% Note: use 'PlaneStrain' for 1D or 2D poroelasticity
Material.constLaw = 'PlaneStrain';

% lumped mass matrix - 0: false, 1: true
Material.lumpedMass = 1;

% lumped damping matrix - 0: false, 1: true
Material.lumpedDamping = 1;

%% Spanos material parameters
% porosity effective pressure coefficient (Spanos, 1989)
% n = 0; % lower limit
n = 1; % return to Biot
% n = Material.Ks/Material.Kf; % upper limit

% modified storage coefficient (Muller, 2019)
Mstarinv = Material.Minv - (1-n)*(Material.alpha - Material.n)/Material.Ks; 
Mstar = 1/Mstarinv;

% porosity equation coefficients
Material.deltaF = (Material.alpha - Material.n) * Material.n * Mstar * n / Material.Ks;
Material.deltaS = (Material.alpha - Material.n) * Material.n * Mstar / Material.Kf;

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
        L = 0.5;
        %%%% solid displacement field
        typeU = 'L3';
        fieldU = 'u';
        MeshU = Build1DMesh(nsd, ne, L, typeU, fieldU);
        %%%% fluid pressure field
        typeP = 'L2';
        fieldP = 'p';
        MeshP = Build1DMesh(nsd, ne, L, typeP, fieldP);
        %%%% porosity field
        if contains(Control.PMmodel, 'UPN')
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
        if contains(Control.PMmodel, 'UPN')
            fieldN = 'n';
            meshFileNameN = 'Column2DQ4.msh';
            MeshN = BuildMesh_GMSH(meshFileNameN, fieldN, nsd, config_dir, progress_on);
        else
            MeshN = [];
        end
end

%% Find nodes for prescribed BCs
% find top and bottom nodes for displacement field
BC.top_node_u = find(MeshU.coords == min(MeshU.coords));
BC.bottom_node_u = find(MeshU.coords == max(MeshU.coords));

% find top and bottom nodes for pressure field
BC.top_node_p = find(MeshP.coords == min(MeshP.coords));
BC.bottom_node_p = find(MeshP.coords == max(MeshP.coords));

%% Initial conditions
% displacement
BC.initU = [];

% pressure
BC.initP = [];

%% Dirichlet BCs - solid
% displacement u=0 at the bottom
BC.fixed_u1 = BC.bottom_node_u;
% u = t at the top
BC.fixed_u2 = BC.top_node_u;
% fixed DOFs
BC.fixed_u = [BC.fixed_u1; BC.fixed_u2];
% fixed DOF values
BC.fixed_u_value1 = zeros(length(BC.fixed_u1),1);
BC.fixed_u_value2 = ones(length(BC.fixed_u2),1);
BC.fixed_u_value = [BC.fixed_u_value1; BC.fixed_u_value2];
% free displacement nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid
%   pressure p=0 at the top
BC.fixed_p = BC.top_node_p;
BC.fixed_p_value = zeros(length(BC.fixed_p),1);
% free pressure nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% point load [GN]
BC.pointLoad = [];

% traction interpolation (needed for traction applied in wells); 1 - true, 0 - false
BC.tractionInterp = 0;

% distributed load [GN/m2]
BC.tractionNodes = [];

% body force [GN/m3]
BC.b = @(x)[];  

%% Neumann BCs - fluid
% point flux [m/s]
BC.pointFlux = [];

% distributed flux [m3/s]
BC.fluxNodes = [];

% flux source [m3/s/m3]
BC.s = @(x)[]; 

%% Porosity BCs
if contains(Control.PMmodel, 'UPN')
    BC.fixed_n = [];
    BC.free_n = setdiff(MeshN.DOF, BC.fixed_n);
    BC.fixed_n_value = zeros(length(BC.fixed_n),1);
end

%% Quadrature order
Control.nqU = 3;
Control.nqP = 3;

%% Frequency domain
Control.freqDomain = 0;  % 1 = true; 0 = false

%% Analytical solution
% 1 = uncoupled problem (elasticity, heat transfer, etc)
% 0 = coupled problem (Biot, Spanos model)
Control.uncoupled = 0; 

% plot analytical solution (valid for 1D problems with Material.Minv == 0)
Control.plotansol = 0; % 1 = true; 0 = false

%% Time step controls
Control.dt = 1.4e-6;  % time step [s]
Control.tend = 1e-3;   % final simulation time [s]

% Newmark method
Control.beta = 0.6;
Control.gamma = 0.6;
Control.theta = 0.6;
Control.lambda = 0.6;

%% Plot data
% DOF to plot graphs
Control.plotu = find(MeshU.coords == 0.1); % x = 0.1m
Control.plotp = find(MeshP.coords == 0.1); % x = 0.1m

end