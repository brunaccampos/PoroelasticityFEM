function [Material, MeshU, MeshP, MeshN, BC, Control] = Plate2D_KirschTest(config_dir, progress_on)
% 2D simulation of hydraulic dilation stimulation of SADG well pair
% Configuration File
% Based on Korsawe (2006) model
% ------------------------------------------------------------------------
% Assumptions/conventions:
% - stress is positive for tension
% - boundary condition for force is based on total stress
% - no acceleration terms for solid or fluid
% - solid velocity is neglected
% - fluid and solid grains are incompressible
% - porosity is constant in space and varies over time
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

% lumped mass matrix - 0: false, 1: true
Material.lumpedMass = 0;

% thickness 
% 1D: cross sectional area [m2]
% 2D: out of plane thickness [m]
Material.t = 1;

% constititive law - 'PlaneStress' or 'PlaneStrain'
% Note: use 'PlaneStrain' for 1D or 2D poroelasticity
Material.constLaw = 'PlaneStress';

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
        ne = 10;
        % column size [m]
        L = 10;
        %%%% solid displacement field
        typeU = 'L3';
        MeshU = Build1DMesh(nsd, ne, L, typeU);
        %%%% fluid pressure field
        typeP = 'L2';
        MeshP = Build1DMesh(nsd, ne, L, typeP);
        %%%% porosity field
        if contains(Control.PMmodel, 'UPN')
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
        meshFileNameU = 'Mesh Files\PlateWithHole_FullGeometryQ4.msh';
        MeshU = BuildMesh_GMSH(meshFileNameU, fieldU, nsd, config_dir, progress_on);
        %%%% pressure field
        fieldP = 'p';
        meshFileNameP = 'Mesh Files\PlateWithHole_FullGeometryQ4.msh';
        MeshP = BuildMesh_GMSH(meshFileNameP, fieldP, nsd, config_dir, progress_on);
        %%%% porosity field
        if contains(Control.PMmodel, 'UPN')
            fieldN = 'n';
            meshFileNameN = 'Mesh Files\PlateWithHole_FullGeometryQ4.msh';
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
% displacement u=0 at boundaries
BC.fixed_u = [MeshU.bottom_dofy; MeshU.top_dofy; MeshU.left_dofx; MeshU.right_dofx];
% fixed DOF values
BC.fixed_u_value = zeros(length(BC.fixed_u),1);
% free nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid
% fixed DOFs (nodes at injection wells)
BC.fixed_p = MeshP.DOF;
% fixed DOF values
BC.fixed_p_value = zeros(length(BC.fixed_p),1);
% free nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% prescribed traction [GN/m2]
BC.tractionNodes = [];

% prescribed traction [GPa]
BC.traction = 10e-3;
% traction interpolation (needed for traction applied in wells); 1 - true, 0 - false
BC.tractionInterp = 1;
% traction unit normal in local coords
normal_L = [0; -1];

% nodes at wells
BC.tractionNodes = [5; 66; 67; 68; 69; 70; 71; 72; 73; 74];

% traction force vector
BC.tractionForceAtNodes = zeros(length(BC.tractionNodes),2);

% coordinates of well centre
centre = [5, 5];

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
    normal_G = rot*normal_L;
    % global traction vector in node i
    BC.tractionForceAtNodes(i,:) = BC.traction * normal_G';
end

% global force vector
BC.tractionForce = zeros(MeshU.nn,2);
for i = 1:length(BC.tractionNodes)
    nodeindex = BC.tractionNodes(i);
    BC.tractionForce(nodeindex,:) = BC.tractionForceAtNodes(i,:);
end

% point loads [GN]
BC.pointLoad = [];

% body force [GN/m3]
BC.b = @(x)[];

%% Neumann BCs - fluid
% distributed flux [m3/s]
% impervious at bottom, left, and right
BC.fluxNodes = [];
BC.fluxValue = zeros(length(BC.fluxNodes),1);

% point flux [m/s]
BC.pointFlux = [];

% flux source [m3/s/m3]
BC.s = @(x)[]; 

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
Control.dt = 1;  % time step [s]
Control.tend = 1;   % final simulation time [s]

Control.beta = 1; % beta-method time discretization -- beta = 1 Backward Euler; beta = 0.5 Crank-Nicolson

%% Plot data
% DOF to plot graphs
Control.plotu = 5*2; % dof y of node 354 (x = 60m, y = 25m)
Control.plotp = 5; % dof of node 209 (x = 60m, y = 25m)

end