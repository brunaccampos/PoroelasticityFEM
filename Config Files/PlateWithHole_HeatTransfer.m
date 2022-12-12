function [Material, MeshU, MeshP, MeshN, BC, Control] = PlateWithHole_HeatTransfer(config_dir, progress_on)
% Plate with hole 1/8 model
% Heat transfer problem adapted from file Q4one8thModel
% ------------------------------------------------------------------------

%% Material properties
% thermal conductance coefficient
Material.kf = 5;

% elasticity modulus [Pa]
Material.E = 0;
% Poisson's ratio
Material.nu = 0;
% Biot's coefficient
Material.alpha = 0;
% 1/Q (related to storage coefficient)
Material.Minv = 0;
% poroelasticity model
Control.Biotmodel = 1;

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
        MeshU.nn = 7; % number of nodes
        MeshU.ne = 3; % number of elements
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
        x1 = -9/4; x2 = -1/4; x3=-1/4*cos(pi/4); x4 = -9/4;
        x5 = 1/2*(x3+x4); x6 =1/2*(x1+x2); x7 = -9/4;
        MeshU.coords(:,1) = [x1;x2;x3;x4;x5;x6;x7];
        % y coordinates
        y1 =0; y2=0; y3 = 1/4*sin(pi/4) ; y4=9/4; 
        y5=1/2*(y3+y4); y6=1/2*(y5+y1); y7 = y6;
        MeshU.coords(:,2) = [y1;y2;y3;y4;y5;y6;y7];

        MeshU.conn = [1,2,6,7; 
                    2,3,5,6;
                    7,6,5,4]; % elements connectivity
                
                
        % mesh for fluid field
        MeshP.nsd = 2; % number of spatial directions
        MeshP.nn = 7; % number of nodes
        MeshP.ne = 3; % number of elements
        MeshP.type = 'Q4'; % element type
        MeshP.field = 'p'; % field type
        MeshP.nne = 4; % nodes per element
        MeshP.nDOFe = MeshP.nne; % DOFs per element
        MeshP.nDOF = MeshP.nn; % total number of DOFs
        MeshP.DOF = (1:MeshP.nDOF).'; % DOFs
        MeshP.coords = zeros(MeshP.nn, 2); % nodal coordinates
        % x coordinates
        x1 = -9/4; x2 = -1/4; x3=-1/4*cos(pi/4); x4 = -9/4;
        x5 = 1/2*(x3+x4); x6 =1/2*(x1+x2); x7 = -9/4;
        MeshP.coords(:,1) = [x1;x2;x3;x4;x5;x6;x7];
        % y coordinates
        y1 =0; y2=0; y3 = 1/4*sin(pi/4) ; y4=9/4; 
        y5=1/2*(y3+y4); y6=1/2*(y5+y1); y7 = y6;
        MeshP.coords(:,2) = [y1;y2;y3;y4;y5;y6;y7];

        MeshP.conn = [1,2,6,7; 
                    2,3,5,6;
                    7,6,5,4]; % elements connectivity

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

%% Dirichlet BCs - solid
% column vector of prescribed displacement dof
BC.fixed_u = 1:MeshU.nDOF;
% prescribed displacement for each dof [u1; u2; ...] [m]
BC.fixed_u_value = zeros(length(BC.fixed_u),1);
% free nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid
BC.fixed_p = [1;2;3;4;7]; % nodes at the inner circle and left boundary
% fixed DOF values
BC.fixed_p_value = [25;10;10;25;25];
% free nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% distributed traction [N/m2]
BC.tractionNodes = [];

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

%% Problem type
% 1 = quasi-steady/transient problem (no acceleration and pressure change)
% 0 = dynamic problem (acceleration/intertia terms included)
Control.steady = 1;

%% Solution parameters
Control.dt = 1;  % time step
Control.tend = 1;   % final simulation time

Control.beta = 1; % beta-method time discretization -- beta = 1 Backward Euler; beta = 0.5 Crank-Nicolson

Control.plotu = 1;
Control.plotp = 1;

% plot analytical solution (valid for 1D problems with Material.Minv == 0)
Control.plotansol = 0; % 1 = true; 0 = false

% solve in the frequency domain
Control.freqDomain = 0;  % 1 = true; 0 = false

end