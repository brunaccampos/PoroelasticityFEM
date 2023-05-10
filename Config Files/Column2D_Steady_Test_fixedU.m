function [Material, MeshU, MeshP, MeshN, BC, Control] = Column2D_Steady_Test_fixedU(config_dir, progress_on)
% Column Consolidation 2D simulation
% Configuration File
% Based on Korsawe (2006) model
% ------------------------------------------------------------------------
% Assumptions/conventions:
% - stress is positive for tension
% - boundary condition for force is based on total stress
% - only solid acceleration is considered (undrained condition; no motions
% of the fluid relative to the solid skeleton can occur)
% - solid grains and fluid are incompressible
% ------------------------------------------------------------------------
% column top at x=L, column bottom at x=0
% ------------------------------------------------------------------------

%% Poroelasticity model
% Options:  Transient_Biot ----- Biot model (u-p), transient
%           Transient_Spanos --- Spanos model (u-p-n), transient
%           Transient_BiotPoro - Biot model (u-p), dynamic, implicit
%                                   porosity perturbation equation
%           Dynamic_Biot ------- Biot model (u-p), dynamic
%           Dynamic_Spanos ----- Spanos model (u-p-n), dynamic
%           Dynamic_BiotPoro --- Biot model (u-p), dynamic, implicit
%                                   porosity perturbation equation
Control.Biotmodel = 'Transient_Biot';

%% Material properties - Komijani (2019)
% shear modulus [GPa]
Material.G = 6000e-3;
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
Material.Kf = 3000e-3;
% solid bulk modulus [GPa]
Material.Ks = 36000e-3;
% fluid bulk viscosity [GPa s]
Material.xif = 2.8e-12; % (Quiroga-Goode, 2005)
% material porosity
Material.n = 0.19;
% 1/Q (related to storage coefficient)
Material.Minv = (Material.alpha - Material.n)/Material.Ks + Material.n/Material.Kf;

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
        ne = 100;
        % column size [m]
        L = 6;
        %%%% solid displacement field
        typeU = 'L3';
        MeshU = Build1DMesh(nsd, ne, L, typeU);
        %%%% fluid pressure field
        typeP = 'L2';
        MeshP = Build1DMesh(nsd, ne, L, typeP);
        %%%% porosity field
        if contains(Control.Biotmodel, 'Spanos')
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
        meshFileNameU = 'Mesh Files\Column2DQ9_Steady_Boone.msh';
        MeshU = BuildMesh_GMSH(meshFileNameU, fieldU, nsd, config_dir, progress_on);
        %%%% pressure field
        fieldP = 'p';
        meshFileNameP = 'Mesh Files\Column2DQ4_Steady_Boone.msh';
        MeshP = BuildMesh_GMSH(meshFileNameP, fieldP, nsd, config_dir, progress_on);
        %%%% porosity field
        if contains(Control.Biotmodel, 'Spanos')
            fieldN = 'n';
            meshFileNameN = 'Mesh Files\Column2DQ4_Steady_Boone.msh';
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

%% Dirichlet BCs - solid
% BC.uy = @(x) x;
% displacement fixed at all nodes
% BC.fixed_u1 = MeshU.DOF(:,1);
% BC.fixed_u2 = MeshU.DOF(:,2);
BC.fixed_u = 1:MeshU.nDOF;
% fixed DOF values
BC.fixed_u_value = zeros(length(BC.fixed_u),1);
% BC.fixed_u_value1 = zeros(length(BC.fixed_u1),1); % zero in x
% BC.fixed_u_value2 = BC.uy(MeshU.coords(:,2)); % linear in y
% 
% BC.fixed_u = [BC.fixed_u1; BC.fixed_u2];
% BC.fixed_u_value = [BC.fixed_u_value1; BC.fixed_u_value2];
% 
% free nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid

% ---------- test 1: all boundary pressures fixed (represent rigid body
% motion)
BC.fixed_p_dof1 = MeshP.left_dof;
BC.fixed_p_dof2 = MeshP.right_dof;
BC.fixed_p_dof3 = MeshP.bottom_dof;
BC.fixed_p_dof4 = MeshP.top_dof;
BC.fixed_p = unique([BC.fixed_p_dof1;BC.fixed_p_dof2;BC.fixed_p_dof3;BC.fixed_p_dof4]);
% fixed DOF values
BC.fixed_p_value = ones(length(BC.fixed_p),1);

% ---------- test 2: linear pressure in x and y (represent constant strain
% state)
% BC.p = @(x) x(:,1) + x(:,2);
% BC.fixed_p = MeshP.DOF;
% BC.fixed_p_value = BC.p(MeshP.coords(BC.fixed_p, :));

% ---------- test 3: linear pressure in x (represent constant strain
% state)
% BC.p = @(x) x(:,1);
% BC.fixed_p = MeshP.DOF;
% BC.fixed_p_value = BC.p(MeshP.coords(BC.fixed_p, 1));

% ---------- test 4: linear pressure in y (represent constant strain
% state)
% BC.p = @(x) x;
% BC.fixed_p = MeshP.DOF;
% BC.fixed_p_value = BC.p(MeshP.coords(BC.fixed_p, 2));

% ---------- test 5: linear pressure in x and y with constant (represent constant strain
% state)
% BC.p = @(x) 1 + 2* x(:,1) + 3*x(:,2);
% BC.fixed_p = MeshP.DOF;
% BC.fixed_p_value = BC.p(MeshP.coords(BC.fixed_p, :));

% ---------- test 6:  pressure constant (represent rigid body motion)
% BC.fixed_p = MeshP.DOF;
% BC.fixed_p_value = ones(length(BC.fixed_p),1);

% ---------- test 7:  pressure zero at bottom
% BC.fixed_p = MeshP.bottom_nodes;
% BC.fixed_p_value = zeros(length(BC.fixed_p),1);

% ---------- test 8:  pressure zero at left
% BC.fixed_p = MeshP.left_nodes;
% BC.fixed_p_value = zeros(length(BC.fixed_p),1);

% ---------- end of testing -----------------------------------------------

% free nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
BC.tractionNodes = [];

% point loads [GN]
BC.pointLoad = [];

% body force [GN/m3]
BC.b = @(x)[];  

%% Neumann BCs - fluid
BC.fluxNodes = [];

% distributed flux [m/s]
% BC.flux = 100;
% 
% % ---------- test 1: flux at top nodes
% BC.fluxNodes = MeshP.top_nodes;
% Flux = BC.flux * max(MeshP.coords(:,1))/length(MeshP.top_nodes);
% BC.fluxValue = Flux * ones(length(BC.fluxNodes),1);
% 
% ---------- test 2: flux at right nodes
% BC.fluxNodes = MeshP.right_nodes;
% Flux = BC.flux * max(MeshP.coords(:,2))/length(MeshP.right_nodes);
% BC.fluxValue = Flux * ones(length(BC.fluxNodes),1);
% ---------- end of testing -----------------------------------------------

% point flux [m3/s]
BC.pointFlux = [];

% flux source [m3/s/m3]
BC.s = @(x)[]; 

%% Quadrature order
Control.nqU = 2;
Control.nqP = 2;

%% Solution parameters
% tag used for computing analytical solution
% 1 = uncoupled problem (elasticity, heat transfer, etc)
% 0 = coupled problem (Biot, Spanos model)
Control.uncoupled = 0; 

Control.dt = 1e-2;  % time step
Control.tend = 10;   % final simulation time

Control.beta = 1; % beta-method time discretization -- beta = 1 Backward Euler; beta = 0.5 Crank-Nicolson

Control.plotu = 1; % dof y of node 63 (x = 0.05m, y = 3m)
Control.plotp = 1; % dof of node 35 (x = 0.05m, y = 3m)

% plot analytical solution (valid for 1D problems with Material.Minv == 0)
Control.plotansol = 0; % 1 = true; 0 = false

% solve in the frequency domain
Control.freqDomain = 0;  % 1 = true; 0 = false

end