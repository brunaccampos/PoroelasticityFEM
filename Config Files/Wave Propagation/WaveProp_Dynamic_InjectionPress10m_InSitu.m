function [Material, MeshU, MeshP, MeshN, BC, Control] = WaveProp_Dynamic_InjectionPress10m_InSitu(config_dir, progress_on,~,~)
% 2D simulation of injection at a well
% Configuration File
% ------------------------------------------------------------------------
% Based on Zienkiewicz (1982) model for dynamic case
% ------------------------------------------------------------------------
% Assumptions/conventions:
% - stress is positive for tension
% - boundary condition for force is based on total stress
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
Control.PMmodel = 'Dyn5_Spanos_UPU';

%% Material properties - Berea Sandstone (Detournay, 1993, p.26)
% elasticity modulus [GPa]
Material.E = 14.4;
% Poisson's ratio
Material.nu = 0.2;
% intrinsic permeability [m2]
Material.k = 1.88e-13;
% dynamic viscosity [GPa s]
Material.mu = 1e-12;
% porous media permeability [m2/GPa s]
Material.kf = Material.k/Material.mu;
% Biot's coefficient
Material.alpha = 0.79;
% fluid bulk modulus [GPa]
Material.Kf = 3.3;
% solid bulk modulus [GPa]
Material.Ks = 36;
% material porosity
Material.n = 0.19;
% 1/Q (related to storage coefficient)
Material.Minv = (Material.alpha - Material.n)/Material.Ks + Material.n/Material.Kf;
% fluid bulk viscosity [GPa s]
Material.xif = 2.8e-12; % (Quiroga-Goode, 2005)
% fluid density [10^9 kg/m3]
Material.rho_f = 1000e-9;
% solid density [10^9 kg/m3]
Material.rho_s = 2600e-9;
% average density of the medium
Material.rho = Material.n*Material.rho_f + (1-Material.n)*Material.rho_s;

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
%%%% pressure field
fieldP = 'p';
meshFileNameP = 'Mesh Files\WavePropInj10x10mQ4_structured.msh';
MeshP = BuildMesh_GMSH(meshFileNameP, fieldP, nsd, config_dir, progress_on);
%%%% porosity field
MeshN = [];

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

% overburden density [10^9 kg/m3]
rho_Top = 2600e-9;
% gravitational acceleration [m/s2]
g = 9.81;
% depth [m]
depth = 1000;
% hydrostatic pressire [GPa]
p0 = Material.rho_f * g * depth;
% vertical in situ stress [GPa]
sigmaV = rho_Top * g * depth;
% horizontal in situ stress [GPa]
sigmaH = sigmaV * Material.nu / (1-Material.nu);
% global in situ stress
sigmaG = [sigmaH, 0; 0, sigmaV];
% traction interpolation (needed for traction applied in wells); 1 - true, 0 - false
BC.tractionInterp = 1;
% traction unit normal in local coords
normalL = [0; 1];

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

% period [s]
t0 = 1e-3;
% peak frequency [Hz]
f = 1/t0;
% amplitude [GN]
a0 = 1e-2;

% adding in situ pressure and injection pressure
BC.tractionForce = @(t) -BC.tractionForce + Material.alpha*p0*normalG_vec + Material.alpha*a0/2*normalG_vec*(1-cos(2*pi*f*t)).*(t<t0);

% point loads [GN]
BC.pointLoad = [];

% body force [GN/m3]
BC.b = @(x,t)[];

%% Neumann BCs - fluid
% point flux [m/s]
BC.pointFlux = [];

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
% 1 = uncoupled problem (elasticity, heat transfer, etc)
% 0 = coupled problem (Biot, Spanos model)
Control.uncoupled = 0; 

% plot analytical solution (valid for 1D problems with Material.Minv == 0)
Control.plotansol = 0; % 1 = true; 0 = false

%% Time step controls
Control.dt = 1e-5;  % time step [s]
Control.tend = 3e-3;   % final simulation time [s]

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