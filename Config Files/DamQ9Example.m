function [Material, MeshU, MeshP, MeshN, BC, Control] = DamQ9Example(config_dir, progress_on)
% Dam example - Q9 element
% Elasticity problem adapted from file Q9_Example_dam_v2
% ------------------------------------------------------------------------

%% Material properties
% elasticity modulus [Pa]
Material.E = 40e9;
% Poisson's ratio
Material.nu = 0.2;
% material density [kg/m3]
Material.rho_c = 2700;
% water density [kg/m3]
Material.rho_w = 1000;

% thermal conductance coefficient
Material.kf = 0;
% Biot's coefficient
Material.alpha = 0;
% 1/Q (related to storage coefficient)
Material.Minv = 0;
% poroelasticity model
Control.Biotmodel = 1;

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
        MeshU.nn = 9; % number of nodes
        MeshU.ne = 1; % number of elements
        MeshU.type = 'Q9'; % element type
        MeshU.field = 'u'; % field type
        MeshU.nne = 9; % nodes per element
        MeshU.nDOFe = MeshU.nne*MeshU.nsd; % DOFs per element
        MeshU.nDOF = MeshU.nn*MeshU.nsd; % total number of DOFs
        for sd = 1:MeshU.nsd
            MeshU.DOF(:,sd) = (sd : MeshU.nsd : (MeshU.nDOF-(MeshU.nsd-sd)))';
        end
        MeshU.coords = zeros(MeshU.nn, 2); % nodal coordinates
        % x coordinates
        MeshU.coords(:,1) = [-1; 1; 0.5; -0.5;0;0.75;0;-0.75;0];
        % y coordinates
        MeshU.coords(:,2) = [0;0;5;5;0;2.5;5;2.5;2.5];
        % elements connectivity
        MeshU.conn = (1:9);       
                
        % mesh for fluid field
        MeshP.nsd = 2; % number of spatial directions
        MeshP.nn = 9; % number of nodes
        MeshP.ne = 1; % number of elements
        MeshP.type = 'Q9'; % element type
        MeshP.field = 'p'; % field type
        MeshP.nne = 9; % nodes per element
        MeshP.nDOFe = MeshP.nne; % DOFs per element
        MeshP.nDOF = MeshP.nn; % total number of DOFs
        MeshP.DOF = (1:MeshP.nDOF).'; % DOFs
        MeshP.coords = zeros(MeshP.nn, 2); % nodal coordinates
        % x coordinates
        MeshP.coords(:,1) = [-1; 1; 0.5; -0.5;0;0.75;0;-0.75;0];
        % y coordinates
        MeshP.coords(:,2) = [0;0;5;5;0;2.5;5;2.5;2.5];
        % elements connectivity
        MeshP.conn = (1:9);       

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
BC.fixed_u = [1;2;3;4;9;10];
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
% distributed traction [N/m2]
BC.tractionNodes = [];

% point loads [N]
BC.pointLoad = zeros(MeshU.nDOF,1);

tractionNodes = [2; 6; 3];
BC.tractionForce = zeros(BC.tractionNodes, 2);
% length where traction is applied
le = sqrt(5^2+0.5^2);
% normal vector at the traction surface
nx = 5/le; % normal x component
ny = 0.5/le; % normal y component
% y coordinates
y = MeshU.coords(tractionNodes,2);
% DOFs traction nodes
dofNodes = [3; 4; 11; 12; 5; 6];
% number of integration points
nq = 3;
% quadrature points
pt = [0.774596669241483; -0.774596669241483; 0.000000000000000];
% quadrature weights
wt = [0.555555555555556; 0.555555555555556; 0.888888888888889];

fe = zeros(6,1);
% loop over integration points
for ip = 1:nq
   csi = pt(ip);
   w = wt(ip);
   % N matrix entries
   N1 = 1/2*(csi)*(csi-1);
   N2 = (1-csi)*(csi+1);
   N3 = 1/2*(csi)*(csi+1);
   N = zeros(2,6);
   N(1,1) = N1;
   N(2,2) = N1;
   N(1,3) = N2;
   N(2,4) = N2;
   N(1,5) = N3;
   N(2,6) = N3;
   % y coordinate at ip
   yip = [N1, N2, N3] * y;
   % traction value at ip
   tbar = - Material.rho_w * 9.81 * (5-yip) * [nx; ny];
   % assemble element force matrix
   fe = fe + N.' * tbar * w * (le/2);
end

BC.pointLoad(dofNodes) = fe;

% body force [N/m3]
BC.b = @(x)[0;
    - Material.rho_c * 9.81];  

%% Neumann BCs - fluid
BC.fluxNodes = [];

% point flux [m/s]
BC.pointFlux = [];

% flux source
BC.s = @(x)[]; 

%% Quadrature order
Control.nqU = 2;
Control.nqP = 2;

%% Problem type
% 1 = quasi-steady/transient problem (no acceleration and pressure change)
% 0 = dynamic problem (acceleration/intertia terms included)
Control.steady = 1;

% tag used for computing analytical solution
% 1 = uncoupled problem (elasticity, heat transfer, etc)
% 0 = coupled problem (Biot, Spanos model)
Control.uncoupled = 0; 

%% Solution parameters
Control.dt = 1;  % time step
Control.tend = 1;   % final simulation time

Control.beta = 1; % beta-method time discretization -- beta = 1 Backward Euler; beta = 0.5 Crank-Nicolson

Control.plotu = 2;
Control.plotp = 2;

% plot analytical solution (valid for 1D problems with Material.Minv == 0)
Control.plotansol = 0; % 1 = true; 0 = false

% solve in the frequency domain
Control.freqDomain = 0;  % 1 = true; 0 = false

end