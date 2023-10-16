function [Material, MeshU, MeshP, MeshN, BC, Control] = WaveProp_Dynamic_InjectionFlux(config_dir, progress_on)
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
Control.PMmodel = 'Dyn4_Biot_UPU';

%% Material properties - Berea Sandstone (Detournay, 1993, p.26)
% elasticity modulus [GPa]
Material.E = 14.4e9;
% Poisson's ratio
Material.nu = 0.2;
% intrinsic permeability [m2]
Material.k = 1.88e-13;
% dynamic viscosity [GPa s]
Material.mu = 1e-3;
% porous media permeability [m2/GPa s]
Material.kf = Material.k/Material.mu;
% Biot's coefficient
Material.alpha = 0.79;
% fluid bulk modulus [GPa]
Material.Kf = 3.3e9;
% solid bulk modulus [GPa]
Material.Ks = 36e9;
% material porosity
Material.n = 0.19;
% 1/Q (related to storage coefficient)
Material.Minv = (Material.alpha - Material.n)/Material.Ks + Material.n/Material.Kf;
% fluid bulk viscosity [GPa s]
Material.xif = 2.8e-3; % (Quiroga-Goode, 2005)
% fluid density [10^9 kg/m3]
Material.rho_f = 1000;
% solid density [10^9 kg/m3]
Material.rho_s = 2600;
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
meshFileNameU = 'Mesh Files\WavePropInj10x10mQ4_finer.msh';
MeshU = BuildMesh_GMSH(meshFileNameU, fieldU, nsd, config_dir, progress_on);
%%%% pressure field
fieldP = 'p';
meshFileNameP = 'Mesh Files\WavePropInj10x10mQ4_finer.msh';
MeshP = BuildMesh_GMSH(meshFileNameP, fieldP, nsd, config_dir, progress_on);
%%%% porosity field
if contains(Control.PMmodel, 'UPN')
    fieldN = 'n';
    meshFileNameN = 'Mesh Files\WavePropInj10x10mQ4_finer.msh';
    MeshN = BuildMesh_GMSH(meshFileNameN, fieldN, nsd, config_dir, progress_on);
else
    MeshN = [];
end

%% Dirichlet BCs - solid
% displacement u=0 normal to all boundaries
BC.fixed_u = [MeshU.left_dofx; MeshU.right_dofx; MeshU.top_dofy; MeshU.bottom_dofy];
% fixed DOF values
BC.fixed_u_value = @(t) zeros(length(BC.fixed_u),1);
% free nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid displacement
% displacement u=0 normal to all boundaries
% BC.fixed_uf = [MeshU.left_dofx; MeshU.right_dofx; MeshU.top_dofy; MeshU.bottom_dofy];
% BC.fixed_uf_value = @(t) zeros(length(BC.fixed_uf),1);

% nodesWell = 4:40;
nodesWell = 4:148; % not transfinite, finer
% peak frequency [Hz]
f = 20;
% peak location [s]
t0 = 1/f;
% amplitude [N]
a0 = 1e3;
BC.fixed_uf = [nodesWell*2-1, nodesWell*2];
BC.fixed_uf_value = @(t) a0*(t-t0).*exp(-(pi*f*(t-t0)).^2).*ones(length(BC.fixed_uf),1);
t = 0:0.001:0.5;
plot(t, BC.fixed_uf_value(t));
% free displacement nodes
BC.free_uf = setdiff(MeshU.DOF, BC.fixed_uf);

%% Dirichlet BCs - fluid
% % nodes at injection well
% % % MeshP.nodesWell = 4:24; % coarse
% % % MeshP.nodesWell = 4:44; % fine
% MeshP.nodesWell = 4:40; % not transfinite
% % % MeshP.nodesWell = 4:76; % not transfinite, fine
% % % MeshP.nodesWell = 4:148; % not transfinite, finer
% % fixed DOFs 
% BC.fixed_p = MeshP.nodesWell;
% % peak frequency [Hz]
% f = 20;
% % peak location [s]
% t0 = 1/f;
% % amplitude [N]
% a0 = 1;
% % fixed DOF values
% BC.fixed_p_value = @(t) a0*(1-2*(pi*f*(t-t0)).^2) .* exp(-(pi*f*(t-t0)).^2)*ones(length(BC.fixed_p),1);

BC.fixed_p = [MeshP.right_nodes; MeshP.top_nodes];
BC.fixed_p_value = @(t) zeros(length(BC.fixed_p),1);
% free nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% point loads [GN]
BC.pointLoad = [];

% prescribed traction [GN/m2]
BC.tractionNodes = [];

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
Control.dt = 1e-4;  % time step [s]
Control.tend = 5e-1;   % final simulation time [s]

% Newmark method
Control.beta = 0.7;
Control.gamma = 0.7;
Control.theta = 0.7;
Control.lambda = 0.7;

%% Plot data
% DOF to plot graphs
% node = 10; % corase
% node = 25; % fine
% node = 14; % not transfinite
% node = 40; % not transfinite, fine
node = 72; % not transfinite, finer

Control.plotu = node*2; % node at the well
Control.plotp = node; % node at the well

end