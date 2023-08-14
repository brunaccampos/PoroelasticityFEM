function [Material, MeshU, MeshP, MeshN, BC, Control] = InjectionWells2D_15x15(config_dir, progress_on)
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
% 1/M (inverse of storage coefficient)
Material.Minv = (Material.alpha - Material.n)/Material.Ks + Material.n/Material.Kf;

% lumped mass matrix - 0: false, 1: true
Material.lumpedMass = 0;

% thickness 
% 1D: cross sectional area [m2]
% 2D: out of plane thickness [m]
Material.t = 1;

% constititive law - 'PlaneStress' or 'PlaneStrain'
% Note: use 'PlaneStrain' for 1D or 2D poroelasticity
Material.constLaw = 'PlaneStrain';

%% Injection data
% water density [kg/m3]
Material.rhof = 1000e-9;
% overburden density [10^9 kg/m3]
Material.rho_Top = 2300e-9;
% gravitational acceleration [m/s2]
Material.g = 9.81;
% depth of the domain [m]
depth = 400;
% hydrostatic pressure [GPa]
p0 = Material.rhof * Material.g * depth;

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
        meshFileNameU = 'Mesh Files\InjectionWells_15x15_Q9.msh';
        MeshU = BuildMesh_GMSH(meshFileNameU, fieldU, nsd, config_dir, progress_on);
        %%%% pressure field
        fieldP = 'p';
        meshFileNameP = 'Mesh Files\InjectionWells_15x15_Q4.msh';
        MeshP = BuildMesh_GMSH(meshFileNameP, fieldP, nsd, config_dir, progress_on);
        %%%% porosity field
        if contains(Control.PMmodel, 'UPN')
            fieldN = 'n';
            meshFileNameN = 'Mesh Files\InjectionWells_15x15_Q4.msh';
            MeshN = BuildMesh_GMSH(meshFileNameN, fieldN, nsd, config_dir, progress_on);
        else
            MeshN = [];
        end
end

%% Initial conditions
% displacement [m]
BC.initU = [];

% pressure [GPa]
% BC.initP = p0 * ones(MeshP.nDOF,1);
BC.initP = [];

%% Dirichlet BCs - solid
% displacement u=0 at left, right, top, and bottom boundaries
BC.fixed_u = [MeshU.bottom_dofy; MeshU.top_dofy; MeshU.left_dofx; MeshU.right_dofx];
% fixed DOF values
BC.fixed_u_value = @(t) zeros(length(BC.fixed_u),1);
% free nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid
% nodes at injection well 2
MeshP.nodesWell2 = [6; 74; 75; 76; 77; 78; 79; 80]; 
% MeshP.nodesWell2 = [6; 146; 147; 148; 149; 150; 151; 152; 153; 154; 155; 156; 157; 158; 159; 160]; % fine mesh
% fixed DOFs 
BC.fixed_p = [MeshP.nodesWell2; MeshP.left_nodes; MeshP.right_nodes];
% BC.fixed_p = [];
% fixed DOF values
BC.fixed_p_value = @(t) zeros(length(BC.fixed_p),1);
% BC.fixed_p_value(1:length(MeshP.nodesWell2)) = -p0;
% free nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% prescribed traction [GN/m2]
BC.tractionNodes = [];

% in situ stress [GPa]
sigmaV = -Material.rho_Top * Material.g * depth + p0;
sigmaH = (-Material.rho_Top * Material.g * depth)*Material.nu/(1-Material.nu) + p0;
% global in situ stress
sigma_G = [sigmaH, 0; 0, sigmaV];

% pressure in well 1 from prescribed flux, applied as a traction
BC.flux = -1e-3;
traction_q = [-BC.flux/Material.kf, 0; 0, -BC.flux/Material.kf];

% traction interpolation (needed for traction applied in wells); 1 - true, 0 - false
BC.tractionInterp = 1;
% traction unit normal in local coords
normal_L = [0; -1];

% nodes at wells
BC.tractionNodes1 = [5; 131; 132; 133; 134; 135; 136; 137; 138; 139; 140; 141; 142; 143; 144; 145];
BC.tractionNodes2 = [6; 146; 147; 148; 149; 150; 151; 152; 153; 154; 155; 156; 157; 158; 159; 160];

% BC.tractionNodes1 = [5; 274; 259; 275; 260; 276; 261; 277; 262; 278; 263; 279; 264; 280; 265; 281; 266; 282; 267; 283; 268; 284; 269; 285; 270; 286; 271; 287; 272; 288; 273; 289]; % fine mesh
% BC.tractionNodes2 = [6; 305; 290; 306; 291; 307; 292; 308; 293; 309; 294; 310; 295; 311; 296; 312; 297; 313; 298; 314; 299; 315; 300; 316; 301; 317; 302; 318; 303; 319; 304; 320]; % fine mesh

BC.tractionNodes = [BC.tractionNodes1; BC.tractionNodes2];

% traction force vector
BC.tractionForce1 = zeros(length(BC.tractionNodes1),2);
BC.tractionForce2 = zeros(length(BC.tractionNodes2),2);

% coordinates of well centre
centre1 = [7.5, 6.5];
centre2 = [7.5, 8.5];

% loop over traction nodes
for i = 1:length(BC.tractionNodes1)
    node = MeshU.coords(BC.tractionNodes1(i),:);
    a = node(1) - centre1(1);
    b = node(2) - centre1(2);
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
    BC.tractionForce1(i,:) = (sigma_G * normal_G)';
%     BC.tractionForce1(i,:) = (sigma_G * normal_G)' + (traction_q * normal_G)';
end

% loop over traction nodes
for i = 1:length(BC.tractionNodes2)
    node = MeshU.coords(BC.tractionNodes2(i),:);
    a = node(1) - centre2(1);
    b = node(2) - centre2(2);
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
    BC.tractionForce2(i,:) = (sigma_G * normal_G)';
end

% weighting for Q9
switch MeshU.type
    case 'Q9'
    BC.tractionForce1(1:8,:) = BC.tractionForce1(1:8,:)*1/3; 
    BC.tractionForce1(9:end,:) = BC.tractionForce1(9:end,:)*2/3; % midside node
    BC.tractionForce2(1:8,:) = BC.tractionForce2(1:8,:)*1/3;
    BC.tractionForce2(9:end,:) = BC.tractionForce2(9:end,:)*2/3; % midside node
end

% traction force vector
BC.tractionForceAtNodes = [BC.tractionForce1; BC.tractionForce2];

% global force vector
BC.tractionForce = zeros(MeshU.nn,2);
for i = 1:length(BC.tractionNodes)
    nodeindex = BC.tractionNodes(i);
    BC.tractionForce(nodeindex,:) = BC.tractionForceAtNodes(i,:);
end

% point loads [GN]
BC.pointLoad = [];

% body force [GN/m3]
BC.b = @(x,t)[];

%% Neumann BCs - fluid
% distributed flux [m3/s]
% impervious at bottom; injection rate at well 2

% flux at well 1
BC.flux = -1e-3;
% nodes well 1
BC.fluxNodes = [5; 67; 68; 69; 70; 71; 72; 73];
% BC.fluxNodes = [5; 131; 132; 133; 134; 135; 136; 137; 138; 139; 140; 141; 142; 143; 144; 145]; % fine mesh
% find surface area
surfs1 = zeros(length(BC.fluxNodes), 2);
L = 0;
count = 1;
for e = 1:MeshP.ne
    conn_e = MeshP.conn(e,:);
    for m = 1:length(BC.fluxNodes)
        nodeM = BC.fluxNodes(m);
        for n = 1:length(BC.fluxNodes)
            nodeN = BC.fluxNodes(n);
            if m ~= n && any(conn_e == nodeM) && any(conn_e == nodeN) && (isequal([nodeM;nodeN], [conn_e(1);conn_e(2)]) || isequal([nodeM;nodeN], [conn_e(2);conn_e(3)]) || isequal([nodeM;nodeN], [conn_e(3);conn_e(4)]) || isequal([nodeM;nodeN], [conn_e(4);conn_e(1)]))
                surfs1(count,1) = nodeM;
                surfs1(count,2) = nodeN;
                L = L + sqrt((MeshP.coords(nodeM,1) - MeshP.coords(nodeN,1))^2 +...
                    (MeshP.coords(nodeM,2) - MeshP.coords(nodeN,2))^2);
                count = count + 1;
            end
        end
    end
end
% total flux
Flux = BC.flux * L/length(BC.fluxNodes);
BC.fluxValue = Flux*ones(size(BC.fluxNodes));

% point flux [m/s]
BC.pointFlux = [];

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
% 1 = uncoupled problem (elasticity, heat transfer, etc)
% 0 = coupled problem (Biot, Spanos model)
Control.uncoupled = 0; 

% plot analytical solution (valid for 1D problems with Material.Minv == 0)
Control.plotansol = 0; % 1 = true; 0 = false

%% Time step controls
Control.dt = 1;  % time step [s]
Control.tend = 50;   % final simulation time [s]

Control.beta = 1; % beta-method time discretization -- beta = 1 Backward Euler; beta = 0.5 Crank-Nicolson

%% Plot data
% DOF to plot graphs
Control.plotu = 1172*2; % point at x=7.46, y=7.45
Control.plotp = 794; % point at x=7.46, y=7.45

% Control.plotu = 1332*2; % fine mesh, point at x=7.49, y=7.53
% Control.plotp = 1172; % fine mesh, point at x=7.49, y=7.53

end