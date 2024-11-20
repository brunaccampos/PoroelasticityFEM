function [Material, MeshU, MeshP, MeshN, BC, Control] = WaveProp_InjPress10m_InSitu(config_dir, progress_on,~,~)
% 2D simulation of injection at a well
% Configuration File
% ------------------------------------------------------------------------
% Based on Zienkiewicz (1982) model for dynamic case
% ------------------------------------------------------------------------
% Assumptions/conventions:
% - stress is positive for tension
% - boundary condition for force is based on total stress
% ------------------------------------------------------------------------
% Porous media theories
% - BT: Biot
% - dCS: de la Cruz and Spanos
% ------------------------------------------------------------------------
% Loading options
% - Tr: transient/quasi-steady
% - Dyn: dynamic (acceleration included)
% ------------------------------------------------------------------------
% Main variables
% u = solid displacement
% p = fluid pressure
% n = porosity
% U = fluid displacement
% v = fluid velocity
% w = relative fluid velocity
% ------------------------------------------------------------------------
% Model options
%
% Tr_BT_UP          Tr_dCS_UP           Tr_dCS_UPN 
%
% Dyn_BT_UP         Dyn_BT_UPU          Dyn_BT_UPV          Dyn_BT_UPW
%
% Dyn_dCS_UP        Dyn_dCS_UPU         Dyn_dCS_UPN         Dyn_dCS_UPW
% ------------------------------------------------------------------------

%% Poroelasticity model
Control.PMmodel = 'Tr_BT_UPU';

%% Material properties - Berea Sandstone (Detournay, 1993, p.26)
% elasticity modulus [GPa]
Material.M(1).E = 14.4;
% Poisson's ratio
Material.M(1).nu = 0.2;
% intrinsic permeability [m2]
Material.M(1).k = 1.88e-13;
% dynamic viscosity [GPa s]
Material.M(1).muf = 1e-12;
% porous media permeability [m2/GPa s]
Material.M(1).kf = Material.M(1).k/Material.M(1).muf;
% Biot's coefficient
Material.M(1).alpha = 0.79;
% fluid bulk modulus [GPa]
Material.M(1).Kf = 3.3;
% solid bulk modulus [GPa]
Material.M(1).Ks = 36;
% material porosity
Material.M(1).eta0 = 0.19;
% 1/Q (related to storage coefficient)
Material.M(1).Minv = (Material.M(1).alpha - Material.M(1).eta0)/Material.M(1).Ks + Material.M(1).eta0/Material.M(1).Kf;
% fluid bulk viscosity [GPa s]
Material.M(1).xif = 2.8e-12; % (Quiroga-Goode, 2005)
% fluid density [10^9 kg/m3]
Material.M(1).rhof = 1000e-9;
% solid density [10^9 kg/m3]
Material.M(1).rhos = 2600e-9;
% average density of the medium
Material.M(1).rho = Material.M(1).eta0*Material.M(1).rhof + (1-Material.M(1).eta0)*Material.M(1).rhos;

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
% n = Material.M(1).Ks/Material.M(1).Kf; % upper limit

% modified storage coefficient (Muller, 2019)
Mstarinv = Material.M(1).Minv - (1-n)*(Material.M(1).alpha - Material.M(1).eta0)/Material.M(1).Ks; 
Mstar = 1/Mstarinv;

% porosity equation coefficients
Material.M(1).deltaf = (Material.M(1).alpha - Material.M(1).eta0) * Material.M(1).eta0 * Mstar * n / Material.M(1).Ks;
Material.M(1).deltas = (Material.M(1).alpha - Material.M(1).eta0) * Material.M(1).eta0 * Mstar / Material.M(1).Kf;

%% Mesh parameters
if progress_on
    disp([num2str(toc),': Building Mesh...']);
end

% Version 2 ASCII - GMSH File
% number of space dimensions
nsd = 2;
%%%% displacement field
fieldU = 'u';
meshFileNameU = 'Mesh Files\WavePropInj10x10mQ9_structured.msh';
MeshU = BuildMesh_GMSH(meshFileNameU, fieldU, nsd, config_dir, progress_on);
% type of material per element
MeshU.MatList = zeros(MeshU.ne, 1, 'int8');
% assign material type to elements
MeshU.MatList(:) = 1;

%%%% pressure field
fieldP = 'p';
meshFileNameP = 'Mesh Files\WavePropInj10x10mQ4_structured.msh';
MeshP = BuildMesh_GMSH(meshFileNameP, fieldP, nsd, config_dir, progress_on);
% type of material per element
MeshP.MatList = zeros(MeshP.ne, 1, 'int8');
% assign material type to elements
MeshP.MatList(:) = 1;

%%%% porosity field
MeshN = [];

%% In situ data
% overburden density [10^9 kg/m3]
rho_Top = 2600e-9;
% gravitational acceleration [m/s2]
g = 9.81;
% depth [m]
depth = 1000;
% hydrostatic pressire [GPa]
p0 = Material.M(1).rhof * g * depth;
% vertical in situ stress [GPa]
sigmaV = rho_Top * g * depth;
% horizontal in situ stress [GPa]
sigmaH = sigmaV * Material.M(1).nu / (1-Material.M(1).nu);
% global in situ stress
sigmaG = [sigmaH, 0; 0, sigmaV];
% traction interpolation (needed for traction applied in wells); 1 - true, 0 - false
BC.tractionInterp = 1;
% traction unit normal in local coords
normalL = [0; 1];

%% Initial BCs
% BC.initP = ones(MeshP.nDOF,1)*p0;

%% Dirichlet BCs - solid
% displacement u=0 at left, right, top, and bottom boundaries
BC.fixed_u = [MeshU.bottom_dofy; MeshU.top_dofy; MeshU.left_dofx; MeshU.right_dofx];
% fixed DOF values
BC.fixed_u_value = @(t) zeros(length(BC.fixed_u),1);
% free nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid displacement
% displacement prescribed on the left and right
BC.fixed_uf = [MeshU.bottom_dofy; MeshU.top_dofy; MeshU.left_dofx; MeshU.right_dofx];
BC.fixed_uf_value = @(t) zeros(length(BC.fixed_uf),1);
% free displacement nodes
BC.free_uf = setdiff(MeshU.DOF, BC.fixed_uf);

%% Dirichlet BCs - fluid
% pressure prescribed at all boundaries
BC.fixed_p = unique([MeshP.left_nodes; MeshP.right_nodes; MeshP.bottom_nodes; MeshP.top_nodes]);
% fixed DOF values
BC.fixed_p_value = @(t) zeros(length(BC.fixed_p),1);

% % pressure prescribed at all boundaries
% BC.fixed_p1 = unique([MeshP.left_nodes; MeshP.right_nodes; MeshP.bottom_nodes; MeshP.top_nodes]);
% BC.fixed_p2 = [1, 5, 12, 321:358]'; % well nodes
% BC.fixed_p = [BC.fixed_p1; BC.fixed_p2];
% % fixed DOF values
% BC.fixed_p_value = @(t) [zeros(length(BC.fixed_p1),1); ones(length(BC.fixed_p2),1)*p0];

% free nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% nodes at well - mesh transfinite, structured
% corner nodes
corner_nodes = [1, 5, 12, 641:659, 680:698];
% middle nodes (Q9 elments)
middle_nodes = [660:679, 699:718];
% all well nodes
BC.tractionNodes = [corner_nodes, middle_nodes];

% traction force vector
BC.tractionAtNodes = zeros(length(BC.tractionNodes),2);
% global unit normal vector
normalG_AtNodes = zeros(length(BC.tractionNodes),2);

% coordinate well centre
centre = [0,0];
% loop over traction nodes
for i = 1:length(BC.tractionNodes)
    node = MeshU.coords(BC.tractionNodes(i),:);
    a = node(1) - centre(1);
    b = node(2) - centre(2);
    theta1 = atan(-a/b);
    theta2 = atan2(-a,b);
    if (a > 0 && b >= 0) % 1st quadrant
        theta = pi()+ theta2;
    elseif (a <= 0 && b > 0) % 2nd quadrant
        theta = pi()+ theta2;
    elseif (a < 0 && b <= 0) % 3rd quadrant
        theta = pi() + theta2;
    elseif (a >= 0 && b < 0) % 4th quadrant
        theta = theta1;
    end
    % rotation matrices for normal vectors
    rot = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    % global normal vector
    normalG = rot*normalL;
    % store normal
    normalG_AtNodes(i,:) = normalG';
    % global traction vector in node i
    BC.tractionAtNodes(i,:) = (sigmaG * normalG)';
    
end

% weighting for Q9
switch MeshU.type
    case 'Q9'
    BC.tractionAtNodes(1:length(corner_nodes),:) = BC.tractionAtNodes(1:length(corner_nodes),:)*1/3; % corner nodes
    BC.tractionAtNodes(length(corner_nodes)+1:end,:) = BC.tractionAtNodes(length(corner_nodes)+1:end,:)*2/3; % midside nodes
end

% global force vector
BC.tractionForce = zeros(MeshU.nn,2);
% vector of unit normals
normalG_vec = zeros(MeshU.nn,2);
% loop over traction nodes
for i = 1:length(BC.tractionNodes)
    nodeindex = BC.tractionNodes(i);
    BC.tractionForce(nodeindex,:) = BC.tractionAtNodes(i,:);
    normalG_vec(nodeindex,:) = normalG_AtNodes(i,:);
end

% adding in situ pressure and injection pressure
BC.tractionForce = @(t) -BC.tractionForce + Material.M(1).alpha*p0*normalG_vec;

% point loads [GN]
BC.pointLoad = @(t)[];

% body force [GN/m3]
BC.b = @(x,t)[];

%% Neumann BCs - fluid
% point flux [m/s]
BC.pointFlux = @(t)[];

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
Control.nqU = 3;
Control.nqP = 3;

%% Frequency domain
Control.freqDomain = 0;  % 1 = true; 0 = false

%% Analytical solution
% plot analytical solution (valid for 1D problems with Material.Minv == 0)
Control.plotansol = 0; % 1 = true; 0 = false

%% Time step controls
Control.dt = 5e-1;  % time step [s]
Control.tend = 100;   % final simulation time [s]

% Newmark method
Control.beta = 0.7;
Control.gamma = 0.7;
Control.theta = 0.7;
Control.lambda = 0.7;

%% Plot data
% DOF to plot graphs (node at the well)
node = 12; % transfinite, structured

Control.plotu = node*2; % y DOF
Control.plotp = node;

end