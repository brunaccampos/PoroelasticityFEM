function [Material, MeshU, MeshP, MeshN, BC, Control] = InjectionWellPoint_Dynamic_Komijani(config_dir, progress_on)
% 2D simulation of hydraulic dilation stimulation of SADG well pair
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
Control.PMmodel = 'Dyn1_Biot_UP';

%% Material properties - Komijani (2019b)
% shear modulus [GPa]
Material.G = 12;
% Poisson's ratio
Material.nu = 0.15;
% elasticity modulus [GPa]
Material.E = 2 * Material.G * (1 + Material.nu);
% intrinsic permeability [m2]
Material.k = 9.869233e-16;
% dynamic viscosity [GPa s]
Material.mu = 1e-12;
% porous media permeability [m2/GPa s]
Material.kf = Material.k/Material.mu;
% Biot's coefficient
Material.alpha = 0.64;
% 1/Q (related to storage coefficient)
Material.Minv = 1/29.3;

% lumped mass matrix - 0: false, 1: true
Material.lumpedMass = 0;

% thickness 
% 1D: cross sectional area [m2]
% 2D: out of plane thickness [m]
Material.t = 1;

% constititive law - 'PlaneStress' or 'PlaneStrain'
% Note: use 'PlaneStrain' for 1D or 2D poroelasticity
Material.constLaw = 'PlaneStrain';

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
        L = 10;
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
        meshFileNameU = 'Mesh Files\Injection20x20 (2)Q9.msh';
        MeshU = BuildMesh_GMSH(meshFileNameU, fieldU, nsd, config_dir, progress_on);
        %%%% pressure field
        fieldP = 'p';
        meshFileNameP = 'Mesh Files\Injection20x20 (2)Q4.msh';
        MeshP = BuildMesh_GMSH(meshFileNameP, fieldP, nsd, config_dir, progress_on);
        %%%% porosity field
        if contains(Control.PMmodel, 'UPN')
            fieldN = 'n';
            meshFileNameN = 'Mesh Files\Injection20x20 (2)Q4.msh';
            MeshN = BuildMesh_GMSH(meshFileNameN, fieldN, nsd, config_dir, progress_on);
        else
            MeshN = [];
        end
end

%% Initial conditions
% displacement
BC.initU = [];

% pressure
p0 = 0; % [GPa]
BC.initP = p0 * ones(MeshP.nDOF,1);

%% Dirichlet BCs - solid
% displacement u=0 at boundaries
BC.fixed_u = [MeshU.bottom_dofy; MeshU.top_dofy; MeshU.left_dofx; MeshU.right_dofx];
% fixed DOF values
BC.fixed_u_value = zeros(length(BC.fixed_u),1);
% free nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid
% fixed DOFs
BC.fixed_p = [MeshP.top_dof; MeshP.bottom_dof; MeshP.left_dof; MeshP.right_dof];
% fixed DOF values
BC.fixed_p_value = zeros(length(BC.fixed_p),1);
% free nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% prescribed traction [GN/m2]
BC.tractionNodes = [];

% point loads [GN]
BC.pointLoad = zeros(MeshU.nDOF,1);
BC.pointLoadNodes = [9; 10];
BC.pointLoad(BC.pointLoadNodes(1)) = 40e-3; % horizontal in situ stress
BC.pointLoad(BC.pointLoadNodes(2)) = 80e-3; % vertical in situ stress

% body force [GN/m3]
BC.b = @(x)[];  

%% Neumann BCs - fluid
% distributed flux [m3/s]
BC.fluxNodes = [];

% point flux [m/s]
BC.pointFluxValue = -8e-4;
BC.pointFluxNodes = 5;
BC.pointFlux = zeros(MeshP.nDOF,1);
BC.pointFlux(BC.pointFluxNodes) = BC.pointFluxValue;

% flux source [m3/s/m3]
BC.s = @(x)[]; 

%% Quadrature order
Control.nqU = 2;
Control.nqP = 2;

%% Frequency domain
Control.freqDomain = 0;  % 1 = true; 0 = false

%% Analytical solution
% 1 = uncoupled problem (elasticity, heat transfer, etc)
% 0 = coupled problem (Biot, Spanos model)
Control.uncoupled = 0; 

% plot analytical solution (valid for 1D problems with Material.Minv == 0)
Control.plotansol = 0; % 1 = true; 0 = false

%% Time step controls
Control.dt = 1;  % time step [s]
Control.tend = 100;   % final simulation time [s]

% Newmark method
Control.beta = 0.7;
Control.gamma = 0.7;
Control.theta = 0.7;

%% Plot graphs
% DOF to plot graphs
Control.plotu = 5*2; % dof y of node 354 (x = 60m, y = 25m)
Control.plotp = 5; % dof of node 209 (x = 60m, y = 25m)

end