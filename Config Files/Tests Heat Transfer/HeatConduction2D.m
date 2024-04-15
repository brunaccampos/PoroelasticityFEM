function [Material, MeshU, MeshP, MeshN, BC, Control] = HeatConduction2D(~, progress_on,~,~)
% Plate with hole 1/8 model: % Heat transfer problem adapted from file Q4one8thModel
% Configuration file
% ------------------------------------------------------------------------

%% Material properties
% thermal conductance coefficient
Material.kf = 5;
% Poroelasticity model
Control.PMmodel = 'Tr1_Biot_UP';

% elasticity modulus [Pa]
Material.E = 0;
% Poisson's ratio
Material.nu = 0;
% Biot's coefficient
Material.alpha = 0;
% 1/Q (related to storage coefficient)
Material.Minv = 0;

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

% Manual 2D mesh
MeshU.nsd = 2; % number of spatial directions
MeshU.nn = 4; % number of nodes
MeshU.ne = 1; % number of elements
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
MeshU.coords(:,1) = [0; 0; 2; 2];
% y coordinates
MeshU.coords(:,2) = [1; 0; 0.5; 1];

MeshU.conn = [1,2,3,4]; % elements connectivity

% mesh for fluid field
MeshP.nsd = 2; % number of spatial directions
MeshP.nn = 4; % number of nodes
MeshP.ne = 1; % number of elements
MeshP.type = 'Q4'; % element type
MeshP.field = 'p'; % field type
MeshP.nne = 4; % nodes per element
MeshP.nDOFe = MeshP.nne; % DOFs per element
MeshP.nDOF = MeshP.nn; % total number of DOFs
MeshP.DOF = (1:MeshP.nDOF).'; % DOFs
MeshP.coords = zeros(MeshP.nn, 2); % nodal coordinates
% x coordinates
MeshP.coords(:,1) = [0; 0; 2; 2];
% y coordinates
MeshP.coords(:,2) = [1; 0; 0.5; 1];

MeshP.conn = [1,2,3,4]; % elements connectivity

% mesh for porosity field
MeshN = [];

%% Dirichlet BCs - solid
% column vector of prescribed displacement dof
BC.fixed_u = 1:MeshU.nDOF;
% prescribed displacement for each dof [u1; u2; ...] [m]
BC.fixed_u_value = @(t) zeros(length(BC.fixed_u),1);
% free nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid
BC.fixed_p = [1; 2; 3]; 
% fixed DOF values
BC.fixed_p_value = @(t) zeros(length(BC.fixed_p),1);
% free nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% distributed traction [N/m2]
BC.tractionNodes = [];

% point loads [N]
BC.pointLoad = [];

% body force [N/m3]
BC.b = @(x,t)[];  

%% Neumann BCs - fluid
BC.fluxNodes = [1; 4];
% distributed flux [m/s]
BC.flux = 20;
Flux = BC.flux *(MeshP.coords(4,1) - MeshP.coords(1,1))/length(BC.fluxNodes);
BC.fluxValue = Flux*ones(size(BC.fluxNodes));

% point flux [m3/s]
BC.pointFlux = [];

% flux source [m3/s/m3]
BC.s = @(x,t) 6; 

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
Control.dt = 1;  % time step
Control.tend = 1;   % final simulation time

Control.beta = 1; % beta-method time discretization -- beta = 1 Backward Euler; beta = 0.5 Crank-Nicolson

%% Plot data
% DOF to plot graphs
Control.plotu = 1;
Control.plotp = 1;

end