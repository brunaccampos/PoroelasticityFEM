function [Material, MeshU, MeshP, MeshN, BC, Control] = Column1D_Steady_Boone(config_dir, progress_on)
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
% column top at x=0, column bottom at x=L
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

%% Material properties - Boone (1990)
% shear modulus [GPa]
Material.G = 6;
% Poisson's ratio
Material.nu = 0.2;
% elasticity modulus [GPa]
Material.E = 2 * Material.G * (1 + Material.nu);
% porous media permeability [m2/GPa s]
Material.kf = 2e-2;
% dynamic viscosity [GPa s]
Material.mu = 1e-12;
% intrinsic permeability [m2]
Material.k = Material.kf * Material.mu;
% Biot's coefficient
Material.alpha = 1;
% fluid bulk modulus [GPa]
Material.Kf = 3;
% solid bulk modulus [GPa]
Material.Ks = 36;
% fluid bulk viscosity [GPa s]
Material.xif = 2.8e-12; % (Quiroga-Goode, 2005)
% material porosity
Material.n = 0.19;
% 1/Q (related to storage coefficient)
Material.Minv = (Material.alpha - Material.n)/Material.Ks + Material.n/Material.Kf;

% additional coefficients for analytical result
% Lame constant [GPa]
Material.lambda = Material.E * Material.nu/((1+Material.nu)*(1-2*Material.nu));
% gravitational acceleration [m/s2]
Material.g = 9.81;
% fluid density [10^9 kg/m3]
Material.rho_f = 1000e-9;
% hydraulic conductivity [m/s]
Material.kh = Material.kf * Material.rho_f * Material.g;

% thickness 
% 1D: cross sectional area [m2]
% 2D: out of plane thickness [m]
Material.t = 1;

% constititive law - 'PlaneStress' or 'PlaneStrain'
% Note: use 'PlaneStrain' for 1D or 2D poroelasticity
Material.constLaw = 'PlaneStrain';

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

% plot deltaS and deltaF
PlotDelta(Material);

% plot coefficients from dimensionless pressure equation
PlotNDPressureEqCoef(Material);

%% Verifying correspondence Biot/Spanos parameters
% pore scale solid constants
Gs = Material.G/(1-Material.n);
lambdaS = Material.Ks - Gs*2/3;
% averaged material constants
M = 1/Material.Minv;
lambda = (1-Material.alpha)*Material.Ks-2*Material.G/3;
% Biot constants
N_Biot = Material.G;
Q_Biot = Material.n * (Material.alpha - Material.n)*M;
R_Biot = Material.n^2*M;
A_Biot = lambda + Q_Biot^2/R_Biot;
% Spanos constants
N_Spanos = (1-Material.n)*Material.G;
Q_Spanos = Material.Ks*Material.deltaF;
R_Spanos = Material.Kf*(Material.n - Material.deltaF);
A_Spanos = (1-Material.n)*lambdaS -Material.deltaS*Material.Ks;
% check term Q^2/R
resBiot = Q_Biot^2/R_Biot;
resSpanos = Q_Spanos^2/R_Spanos;
% print info
fprintf('Biot Constants: \n N = %.4f \n Q = %.4f \n R = %.4f \n A = %.4f \n Q^2/R = %.4f \n', N_Biot, Q_Biot, R_Biot, A_Biot, resBiot);
fprintf('Spanos Constants: \n N = %.4f \n Q = %.4f \n R = %.4f \n A = %.4f \n Q^2/R = %.4f \n', N_Spanos, Q_Spanos, R_Spanos, A_Spanos, resSpanos);

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
        ne = 16;
        % column size [m]
        L = 6;
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
BC.fixed_u = (BC.bottom_node_u);
BC.fixed_u_value = @(t) zeros(length(BC.fixed_u),1);
% free displacement nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid
%   pressure p=0 at the top
BC.fixed_p = (BC.top_node_p);
BC.fixed_p_value = @(t) zeros(length(BC.fixed_p),1);
% free pressure nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% point load [GN]
BC.pointLoadValue = 1e-6;
BC.pointLoadNodes = BC.top_node_u;
BC.pointLoad = zeros(MeshU.nDOF,1);
BC.pointLoad(BC.pointLoadNodes) = BC.pointLoadValue;

% traction interpolation (needed for traction applied in wells); 1 - true, 0 - false
BC.tractionInterp = 0;

% distributed load [GN/m2]
BC.tractionNodes = [];

% body force [GN/m3]
BC.b = @(x,t)[];  

%% Neumann BCs - fluid
% point flux [m/s]
BC.pointFluxValue = 0;
BC.pointFluxNodes = [];
BC.pointFlux = zeros(MeshP.nDOF,1);
BC.pointFlux(BC.pointFluxNodes) = BC.pointFluxValue;

% distributed flux [m3/s]
BC.fluxNodes = [];

% flux source [m3/s/m3]
BC.s = @(x,t)[]; 

%% Porosity BCs
if contains(Control.PMmodel, 'UPN')
    BC.fixed_n = [];
    BC.free_n = setdiff(MeshN.DOF, BC.fixed_n);
    BC.fixed_n_value = zeros(length(BC.fixed_n),1);
end

%% Quadrature order
Control.nqU = 2;
Control.nqP = 1;

%% Frequency domain
Control.freqDomain = 0;  % 1 = true; 0 = false

%% Solution parameters
% tag used for computing analytical solution
% 1 = uncoupled problem (elasticity, heat transfer, etc)
% 0 = coupled problem (Biot, Spanos model)
Control.uncoupled = 0; 

% plot analytical solution (valid for 1D problems with Material.Minv == 0)
Control.plotansol = 1; % 1 = true; 0 = false

%% Time step controls
Control.dt = 1e-4;  % time step [s]
Control.tend = 1;   % final simulation time [s]

Control.beta = 1; % beta-method time discretization -- beta = 1 Backward Euler; beta = 0.5 Crank-Nicolson

%% Plot data
% DOF to plot graphs
Control.plotu = round(length(MeshU.coords)/2);
Control.plotp = round(length(MeshP.coords)/2);

end