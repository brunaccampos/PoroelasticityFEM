function [Material, MeshU, MeshP, MeshN, BC, Control] = Beam_Dynamic(config_dir, progress_on)
% Plate with hole 1/8 model
% Heat transfer problem adapted from file Q4one8thModel
% ------------------------------------------------------------------------

%% Material properties
% shear modulus [Pa]
Material.G = 10;
% Poisson's ratio
Material.nu = 0.3;
% elasticity modulus [Pa]
Material.E = 2 * Material.G * (1 + Material.nu);
% material density [kg/m3]
Material.rho = 1;

% thermal conductance coefficient
Material.kf = 0;
% Biot's coefficient
Material.alpha = 0;
% 1/Q (related to storage coefficient)
Material.Minv = 0;
% fluid density [kg/m3]
Material.rho_f = 0;
% poroelasticity model
Control.Biotmodel = 1;

% constititive law - 'PlaneStress' or 'PlaneStrain'
% Note: use 'PlaneStrain' for 1D or 2D poroelasticity
Material.constLaw = 'PlaneStress';

% lumped mass matrix - 0: false, 1: true
Material.lumpedMass = 0;

%% Mesh Properties
if progress_on
    disp([num2str(toc),': Building Mesh...']);
end

% mesh type
% 'Manual': 1D mesh
% 'Gmsh': 2D mesh, input file from GMSH
MeshType = 'Manual';

switch MeshType
    case 'Manual'
%         % number of space dimensions
%         nsd = 1;
%         % number of elements
%         ne = 10;
%         % column size [m]
%         L = 10;
%         %%%% solid displacement field
%         typeU = 'L3';
%         MeshU = Build1DMesh(nsd, ne, L, typeU);
%         %%%% fluid pressure field
%         typeP = 'L2';
%         MeshP = Build1DMesh(nsd, ne, L, typeP);
%         %%%% porosity field
%         if ~Control.Biotmodel
%             typeN = 'L2';
%             MeshN = Build1DMesh(nsd, ne, L, typeN);
%         else
%             MeshN = [];
%         end
% ------------------------------------------------------------------------
% Manual 2D mesh
        MeshU.nsd = 2; % number of spatial directions
        MeshU.nn = 22; % number of nodes
        MeshU.ne = 10; % number of elements
        MeshU.type = 'Q4'; % element type
        MeshU.field = 'u'; % field type
        MeshU.nne = 4; % nodes per element
        MeshU.nDOFe = MeshU.nne*MeshU.nsd; % DOFs per element
        MeshU.nDOF = MeshU.nn*MeshU.nsd; % total number of DOFs
        for sd = 1:MeshU.nsd
            MeshU.DOF(:,sd) = (sd : MeshU.nsd : (MeshU.nDOF-(MeshU.nsd-sd)))';
        end
        MeshU.coords = zeros(MeshU.nn, 2); % nodal coordinates
        % x coordinates
        MeshU.coords(:,1) = [(0:10)'; (0:10)'];
        % y coordinates
        MeshU.coords(:,2) = [zeros(11,1); ones(11,1)];
        % elements connectivity
        MeshU.conn = [1, 2, 13, 12;
                    2, 3, 14, 13;
                    3, 4, 15, 14;
                    4, 5, 16, 15;
                    5, 6, 17, 16;
                    6, 7, 18, 17;
                    7, 8, 19, 18;
                    8, 9, 20, 19;
                    9, 10, 21, 20;
                    10, 11, 22, 21]; 
                
        % mesh for fluid field
        MeshP.nsd = 2; % number of spatial directions
        MeshP.nn = 22; % number of nodes
        MeshP.ne = 10; % number of elements
        MeshP.type = 'Q4'; % element type
        MeshP.field = 'p'; % field type
        MeshP.nne = 4; % nodes per element
        MeshP.nDOFe = MeshP.nne; % DOFs per element
        MeshP.nDOF = MeshP.nn; % total number of DOFs
        MeshP.DOF = (1:MeshP.nDOF).'; % DOFs
        MeshP.coords = zeros(MeshP.nn, 2); % nodal coordinates
        % x coordinates
        MeshP.coords(:,1) = [(0:10)'; (0:10)'];
        % y coordinates
        MeshP.coords(:,2) = [zeros(11,1); ones(11,1)];
        % elements connectivity
        MeshP.conn = [1, 2, 13, 12;
                    2, 3, 14, 13;
                    3, 4, 15, 14;
                    4, 5, 16, 15;
                    5, 6, 17, 16;
                    6, 7, 18, 17;
                    7, 8, 19, 18;
                    8, 9, 20, 19;
                    9, 10, 21, 20;
                    10, 11, 22, 21]; 

        % mesh for porosity field
        MeshN = [];
% ------------------------------------------------------------------------
    case 'Gmsh'
        % Version 2 ASCII
        % number of space dimensions
        nsd = 2;
        %%%% displacement field
        fieldU = 'u';
        meshFileNameU = 'PlateWithHoleQ4.msh';
        MeshU = BuildMesh_GMSH(meshFileNameU, fieldU, nsd, config_dir, progress_on);
        %%%% pressure field
        fieldP = 'p';
        meshFileNameP = 'PlateWithHoleQ4.msh';
        MeshP = BuildMesh_GMSH(meshFileNameP, fieldP, nsd, config_dir, progress_on);
        %%%% porosity field
        if ~Control.Biotmodel
            fieldN = 'n';
            meshFileNameN = 'PlateWithHoleQ4.msh';
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
% column vector of prescribed displacement dof
BC.fixed_u = [1; 2; 23];
% prescribed displacement for each dof [u1; u2; ...] [m]
BC.fixed_u_value = zeros(length(BC.fixed_u),1);
% free nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid
BC.fixed_p = 1:MeshP.nDOF; % nodes at the inner circle
% fixed DOF values
BC.fixed_p_value = zeros(length(BC.fixed_p),1);
% free nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% traction interpolation (needed for traction applied in wells); 1 - true, 0 - false
BC.tractionInterp = 0;

BC.tractionNodes = [11; 22];
% prescribed traction [N/m2]
BC.traction = 0.005;

Fright = BC.traction*max(MeshU.coords(:,2))/((length(BC.tractionNodes) - 1)/2);
BC.tractionForce = zeros(length(BC.tractionNodes),2);
BC.tractionForce = [zeros(size(BC.tractionNodes)), Fright/2 *ones(size(BC.tractionNodes))];

% find the nodes in the top left and bottom right corners
botrightnode = find(MeshU.coords(BC.tractionNodes,2) == min(MeshU.coords(:,2)));
toprightnode  = find(MeshU.coords(BC.tractionNodes,2) == max(MeshU.coords(:,2)));

BC.tractionForce(botrightnode,2) = BC.tractionForce(botrightnode,2)/2;
BC.tractionForce(toprightnode,2) = BC.tractionForce(toprightnode,2)/2;

% point loads [N]
BC.pointLoad = [];

% body force [N/m3]
BC.b = @(x)[];  

%% Neumann BCs - fluid
% distributed flux [m3/s]
BC.fluxNodes = 1:MeshP.nDOF;
BC.fluxValue = zeros(length(BC.fluxNodes),1);

% point flux [m/s]
BC.pointFlux = [];

% flux source [m3/s/m3]
BC.s = @(x)[]; 

%% Quadrature order
Control.nqU = 2;
Control.nqP = 2;

%% Problem type
% 1 = quasi-steady/transient problem (no acceleration and pressure change)
% 0 = dynamic problem (acceleration/intertia terms included)
Control.steady = 0;

%% Solution parameters
Control.dt = 0.1;  % time step
Control.tend = 500;   % final simulation time
Control.tol = 1e-3; % tolerance for NR method
Control.max_it = 100; % maximum of iterations

Control.plotu = 44;
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