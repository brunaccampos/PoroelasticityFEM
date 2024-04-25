function [Material, MeshU, MeshP, MeshN, BC, Control] = MasterConfigFile(config_dir, progress_on,~,~)
% [PROBLEM TYPE, e.g., column consolidation]
% Configuration File

% [ASSUMPTIONS/CONVENTIONS ACCORDING TO TRANSIENT/DYNAMIC CASE]
% ------------------------------------------------------------------------
% Based on Korsawe (2006) model for transient/quasi-steady case
% ------------------------------------------------------------------------
% Assumptions/conventions:
% - stress is positive for tension
% - boundary condition for force is based on total stress
% - no acceleration terms for solid or fluid
% - solid velocity is neglected
% ------------------------------------------------------------------------

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
Control.PMmodel = 'Tr1_Biot_UP';

%% Material properties - [INSERT REFERENCE]
% Poisson's ratio
Material.nu = 0.2;
% dynamic viscosity [GPa s]
Material.muf = 1e-12;
% intrinsic permeability [m2]
Material.k = 1e-13;
% porous media permeability [m2/GPa s]
Material.kf = Material.k / Material.muf;
% Biot's coefficient
Material.alpha = 1;
% fluid bulk modulus [GPa]
Material.Kf = 2.2;
% solid bulk modulus [GPa]
Material.Ks = 33;
% fluid bulk viscosity [GPa s]
Material.xif = 2.8e-12;
% material porosity
Material.eta0 = 0.25;
% shear modulus [GPa]
Material.mu = (1-Material.eta0)*23;
% elasticity modulus [GPa]
Material.E = 2 * Material.mu * (1 + Material.nu);
% 1/Q (related to storage coefficient)
Material.Minv = (Material.alpha - Material.eta0)/Material.Ks + Material.eta0/Material.Kf;
% fluid density [10^9 kg/m3]
Material.rhof = 1000e-9;
% solid density [10^9 kg/m3]
Material.rhos = 2650e-9;
% average density of the medium
Material.rho = Material.eta0*Material.rhof + (1-Material.eta0)*Material.rhos;
% added mass [10^9 kg/m3]
Material.rho12 = -83e-9;

% thickness 
% 1D: cross sectional area [m2]
% 2D: out of plane thickness [m]
Material.t = 1;

% constititive law - 'PlaneStress' or 'PlaneStrain'
% Note: use 'PlaneStrain' for 1D or 2D poroelasticity
Material.constLaw = 'PlaneStrain';

% lumped mass matrix - 0: false, 1: true
Material.lumpedMass = 0;

% lumped damping matrix - 0: false, 1: true
Material.lumpedDamping = 0;

%% Spanos material parameters
% porosity effective pressure coefficient (Spanos, 1989)
% n = 0; % lower limit
n = 1; % return to Biot
% n = Material.Ks/Material.Kf; % upper limit

% modified storage coefficient (Muller, 2019)
Mstarinv = Material.Minv - (1-n)*(Material.alpha - Material.eta0)/Material.Ks; 
Mstar = 1/Mstarinv;

Material.deltaf = (Material.alpha - Material.eta0) * Material.eta0 * Mstar * n / Material.Ks;
Material.deltas = (Material.alpha - Material.eta0) * Material.eta0 * Mstar /Material.Kf;

% plot deltaS and deltaF
PlotDelta(Material);

% plot coefficients from dimensionless pressure equation
PlotNDPressureEqCoef(Material);

%% Verifying correspondence Biot/Spanos parameters
% pore scale solid constants
mus = Material.mu/(1-Material.eta0);
lambdaS = Material.Ks - mus*2/3;
% averaged material constants
M = 1/Material.Minv;
lambda = (1-Material.alpha)*Material.Ks-2*Material.mu/3;
% Biot constants
N_Biot = Material.mu;
Q_Biot = Material.eta0 * (Material.alpha - Material.eta0)*M;
R_Biot = Material.eta0^2*M;
A_Biot = lambda + Q_Biot^2/R_Biot;
% Spanos constants
N_Spanos = (1-Material.eta0)*Material.mu;
Q_Spanos = Material.Ks*Material.deltaf;
R_Spanos = Material.Kf*(Material.eta0 - Material.deltaf);
A_Spanos = (1-Material.eta0)*lambdaS -Material.deltas*Material.Ks;
% check term Q^2/R
resBiot = Q_Biot^2/R_Biot;
resSpanos = Q_Spanos^2/R_Spanos;
% print info
fprintf('Biot Constants: \n N = %.4f \n Q = %.4f \n R = %.4f \n A = %.4f \n Q^2/R = %.4f \n', N_Biot, Q_Biot, R_Biot, A_Biot, resBiot);
fprintf('Spanos Constants: \n N = %.4f \n Q = %.4f \n R = %.4f \n A = %.4f \n Q^2/R = %.4f \n', N_Spanos, Q_Spanos, R_Spanos, A_Spanos, resSpanos);

%% Mesh parameters
if progress_on
    disp([num2str(toc),': Building Mesh...']);
end

% mesh type: 'Manual' or 'Gmsh'
MeshType = 'Gmsh';

switch MeshType
    case 'Manual'
        % location of initial node [m] [x0;y0;z0]
        coord0 = [0;0;0];
        % number of space dimensions
        nsd = 1;
        % size of domain [m] [Lx;Ly;Lz]
        L = 10;
        % number of elements in each direction [nex; ney; nez]
        ne = 100;
        
        %%%% displacement mesh
        % element type ('Q4')
        typeU = 'L3';
        % variable field ('u', 'p', 'n')
        fieldU = 'u';
        MeshU = BuildMesh_structured(nsd, coord0, L, ne, typeU, fieldU, progress_on);
        
        %%%% pressure mesh
        % element type ('Q4')
        typeP = 'L2';
        % variable field ('u', 'p', 'n')
        fieldP = 'p';
        MeshP = BuildMesh_structured(nsd, coord0, L, ne, typeP, fieldP, progress_on);
        
        %%%% porosity mesh
        if contains(Control.PMmodel, 'UPN')
            % element type ('Q4')
            typeN = 'L2';
            % variable field ('u', 'p', 'n')
            fieldN = 'n';
            MeshN = BuildMesh_structured(nsd, coord0, L, ne, typeN, fieldN, progress_on);
        else
            MeshN = [];
        end
    case 'Gmsh'
        % Version 2 ASCII
        % number of space dimensions
        nsd = 2;
        %%%% displacement field
        fieldU = 'u';
        meshFileNameU = 'Mesh Files\Plate_15x15Q4finer.msh';
        MeshU = BuildMesh_GMSH(meshFileNameU, fieldU, nsd, config_dir, progress_on);
        %%%% pressure field
        fieldP = 'p';
        meshFileNameP = 'Mesh Files\Plate_15x15Q4finer.msh';
        MeshP = BuildMesh_GMSH(meshFileNameP, fieldP, nsd, config_dir, progress_on);
        %%%% porosity field
        if contains(Control.PMmodel, 'UPN')
            fieldN = 'n';
            meshFileNameN = 'Mesh Files\Plate_15x15Q4finer.msh';
            MeshN = BuildMesh_GMSH(meshFileNameN, fieldN, nsd, config_dir, progress_on);
        else
            MeshN = [];
        end
end

%% Dirichlet BCs - solid
% central node
node = find(MeshU.coords(:,1) == 7.5 & MeshU.coords(:,2) == 7.5);
% central node y DOF
BC.fixed_u = node*2;
% period [ms]
t0 = 1e-3;
% fixed DOF values
BC.fixed_u_value = @(t) (sin(2*pi*(t)/t0) - 0.5*sin(4*pi*(t)/t0)).*(t<t0);
% BC.fixed_u_value = @(t) (-(t0/2/pi)*cos(2*pi*(t)/t0) + (t0/8/pi)*cos(4*pi*(t)/t0)).*(t<t0);
% free displacement nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid
BC.fixed_p = [];
% fixed DOF values
BC.fixed_p_value = @(t) zeros(length(BC.fixed_p),1);
% free pressure nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% point load [GN]
BC.pointLoad = @(t)[];

% distributed load [GN/m2]
BC.tractionNodes = [];

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
Control.nqU = 2;
Control.nqP = 2;

%% Frequency domain
Control.freqDomain = 0;  % 1 = true; 0 = false

%% Analytical solution
% plot analytical solution (valid for 1D problems)
Control.plotansol = 0; % 1 = true; 0 = false

% type of analytical solution to compute
% 'getAnSol_uncoupled' = uncoupled problem (elasticity, heat transfer, etc)
% 'getAnSol_coupledComp' = coupled porous media problem, compressible
% materials
% 'getAnSol_coupledIncomp' = coupled porous media problem, incompressible
% materials (1/M=0)
Control.ansol_type = 'getAnSol_uncoupled';

%% Time step controls
Control.dt = 1e-5;  % time step
Control.tend = 2e-3;   % final simulation time

% Newmark method
Control.beta = 0.7;
Control.gamma = 0.7;
Control.theta = 0.7;
Control.lambda = 0.7;

% HHT method
Control.alpha = 0;

% adaptive time step (optional)
% Control.dtmin = 1e-4; % minimum time step
% Control.tlim = 1; % limit to use dtmin

% ramp load option (optional); uses tlim from adaptive time step
% Control.rampLoad = 1;

%% Plot data
% DOF to plot graphs
Control.plotu = node*2; % dof y of node 242 (x = 7.5m, y = 7.5m)
Control.plotp = node; % dof y of node 177 (x = 7.5m, y = 7.5m)

% Plot synthetics
Control.plotSyntheticsON = 1; % 0: false, 1: true

% Plot in a row
Control.fixedDepthPlotON = 1; % 0: false, 1: true

Control.depthplot = 7.5; % fixed coordinate
Control.depthDir = 1; % 1 = fixed y, vary x --- 2 = fixed x, vary y

% node numbering
switch Control.depthDir
    case 1
        rowofnodes_u = find(MeshU.coords(:,2) == Control.depthplot);
        rowofnodes_p = find(MeshP.coords(:,2) == Control.depthplot); 
    case 2
        rowofnodes_u = find(MeshU.coords(:,1) == Control.depthplot); 
        rowofnodes_p = find(MeshP.coords(:,1) == Control.depthplot); 
end

nodes_u = [MeshU.coords(rowofnodes_u,Control.depthDir), rowofnodes_u]; % matrix with node numbering and variable coord
nodes_u_sorted = sortrows(nodes_u); % order in terms of variable coord
Control.ploturow = [nodes_u_sorted(:,2) .* 2 - 1; nodes_u_sorted(:,2) .* 2];

nodes_p = [MeshP.coords(rowofnodes_p,Control.depthDir), rowofnodes_p]; % matrix with node numbering and variable coord
nodes_p_sorted = sortrows(nodes_p); % order in terms of variable coord
Control.plotprow = nodes_p_sorted(:,2);

end