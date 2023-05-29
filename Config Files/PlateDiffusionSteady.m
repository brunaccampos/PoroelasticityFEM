function [Material, MeshU, MeshP, MeshN, BC, Control] = PlateDiffusionSteady(config_dir, progress_on)
% ------------------------------------------------------------------------
% 2D diffusion problem - steady case
% Adapted from: https://github.com/GCMLab
% ------------------------------------------------------------------------

%% Material properties
% diffusion coefficient [m2/s] 
Material.kf = 4e-6;
% 1/Q (related to storage coefficient)
Material.Minv = 0;

% material density [kg/m3]
Material.rho = 0;
% fluid density [kg/m3]
Material.rho_f = 0;
% elasticity modulus [Pa]
Material.E = 0;
% Poisson's ratio
Material.nu = 0;
% Biot's coefficient
Material.alpha = 0;
% poroelasticity model
Control.PMmodel = 'Tr1_Biot_UP';

% lumped mass matrix - 0: false, 1: true
Material.lumpedMass = 0;

% thickness 
% 1D: cross sectional area [m2]
% 2D: out of plane thickness [m]
Material.t = 1;

% constititive law - 'PlaneStress' or 'PlaneStrain'
% Note: use 'PlaneStrain' for 1D or 2D poroelasticity
Material.constLaw = 'PlaneStress';

%% Mesh Properties
if progress_on
    disp([num2str(toc),': Building Mesh...']);
end

% mesh type
% 'Manual': 1D mesh
% 'Gmsh': 2D mesh, input file from GMSH
MeshType = 'Gmsh';

switch MeshType
    case 'Manual'
    % Manual 2D mesh
        MeshU.nsd = 2; % number of spatial directions
        MeshU.nn = 9; % number of nodes
        MeshU.ne = 8; % number of elements
        MeshU.type = 'T3'; % element type
        MeshU.field = 'u'; % field type
        MeshU.nne = 3; % nodes per element
        MeshU.nDOFe = MeshU.nne*MeshU.nsd; % DOFs per element
        MeshU.nDOF = MeshU.nn*MeshU.nsd; % total number of DOFs
        for sd = 1:MeshU.nsd
            MeshU.DOF(:,sd) = (sd : MeshU.nsd : (MeshU.nDOF-(MeshU.nsd-sd)))';
        end
        MeshU.coords = zeros(MeshU.nn, 2); % nodal coordinates
        % x coordinates
        MeshU.coords(:,1) = [0; 0.5; 1; 0; 0.5; 1; 0; 0.5; 1];
        % y coordinates
        MeshU.coords(:,2) = [0; 0; 0; 0.5; 0.5; 0.5; 1; 1; 1];
        MeshU.conn = [1, 5, 4;
            1, 2, 5;
            2, 6, 5;
            2, 3, 6;
            4, 8, 7;
            4, 5, 8;
            5, 9, 8;
            5, 6, 9]; % elements connectivity
                
                
        % mesh for fluid field
        MeshP.nsd = 2; % number of spatial directions
        MeshP.nn = 9; % number of nodes
        MeshP.ne = 8; % number of elements
        MeshP.type = 'T3'; % element type
        MeshP.field = 'p'; % field type
        MeshP.nne = 3; % nodes per element
        MeshP.nDOFe = MeshP.nne; % DOFs per element
        MeshP.nDOF = MeshP.nn; % total number of DOFs
        MeshP.DOF = (1:MeshP.nDOF).'; % DOFs
        MeshP.coords = zeros(MeshP.nn, 2); % nodal coordinates
        % x coordinates
        MeshP.coords(:,1) = [0; 0.5; 1; 0; 0.5; 1; 0; 0.5; 1];
        % y coordinates
        MeshP.coords(:,2) = [0; 0; 0; 0.5; 0.5; 0.5; 1; 1; 1];
        MeshP.conn = [1, 5, 4;
            1, 2, 5;
            2, 6, 5;
            2, 3, 6;
            4, 8, 7;
            4, 5, 8;
            5, 9, 8;
            5, 6, 9]; % elements connectivity

        % mesh for porosity field
        MeshN = [];
% ------------------------------------------------------------------------
    case 'Gmsh'
        % Version 2 ASCII
        % number of space dimensions
        nsd = 2;
        %%%% displacement field
        fieldU = 'u';
        meshFileNameU = 'Mesh Files\PlateDiffusion.msh';
        MeshU = BuildMesh_GMSH(meshFileNameU, fieldU, nsd, config_dir, progress_on);
        %%%% pressure field
        fieldP = 'p';
        meshFileNameP = 'Mesh Files\PlateDiffusion.msh';
        MeshP = BuildMesh_GMSH(meshFileNameP, fieldP, nsd, config_dir, progress_on);
        %%%% porosity field
        if contains(Control.PMmodel, 'UPN')
            fieldN = 'n';
            meshFileNameN = 'Mesh Files\PlateDiffusion.msh';
            MeshN = BuildMesh_GMSH(meshFileNameN, fieldN, nsd, config_dir, progress_on);
        else
            MeshN = [];
        end
end

%% Dirichlet BCs - solid
% column vector of prescribed displacement dof
BC.fixed_u = 1:MeshU.nDOF;
% prescribed displacement for each dof [u1; u2; ...] [m]
BC.fixed_u_value = zeros(length(BC.fixed_u),1);
% free nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid
BC.fixed_p = unique([MeshP.bottom_nodes; MeshP.top_nodes; MeshP.right_nodes; MeshP.left_nodes]);
% fixed DOF values
BC.fixed_p_value = zeros(length(BC.fixed_p),1);
% free nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% distributed traction [N/m2]
BC.tractionNodes = 1:MeshU.nDOF;
BC.tractionForce = zeros(length(BC.tractionNodes),2);

% point loads [N]
BC.pointLoad = [];

% body force [N/m3]
BC.b = @(x)[];  

%% Neumann BCs - fluid
% distributed flux [m/s]
BC.fluxNodes = [];

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

% basic time step controls
Control.dt = 1;  % time step
Control.tend = 1;   % final simulation time

Control.beta = 1; % beta-method time discretization -- beta = 1 Backward Euler; beta = 0.5 Crank-Nicolson

% DOF to plot graphs
Control.plotu = 2;
Control.plotp = 2;

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