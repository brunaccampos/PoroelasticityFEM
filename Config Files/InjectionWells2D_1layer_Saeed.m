function [Material, MeshU, MeshP, MeshN, BC, Control] = InjectionWells2D_1layer_Saeed(config_dir, progress_on)
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
% 1 - Biot theory
% 0 - Spanos theory (additional porosity equation)
Control.Biotmodel = 0;

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
Material.constLaw = 'PlaneStrain';

%% Spanos material parameters
% porosity effective pressure coefficient (Spanos, 1989)
% n = 0; % lower limit
% n = 1; % return to Biot
n = Material.Ks/Material.Kf; % upper limit

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
        if ~Control.Biotmodel
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
        meshFileNameU = 'Mesh Files\InjectionWellsQ9.msh';
        MeshU = BuildMesh_GMSH(meshFileNameU, fieldU, nsd, config_dir, progress_on);
        %%%% pressure field
        fieldP = 'p';
        meshFileNameP = 'Mesh Files\InjectionWellsQ4.msh';
        MeshP = BuildMesh_GMSH(meshFileNameP, fieldP, nsd, config_dir, progress_on);
        %%%% porosity field
        if ~Control.Biotmodel
            fieldN = 'n';
            meshFileNameN = 'Mesh Files\InjectionWellsQ4.msh';
            MeshN = BuildMesh_GMSH(meshFileNameN, fieldN, nsd, config_dir, progress_on);
        else
            MeshN = [];
        end
end

%% Initial conditions
% displacement
BC.initU = [];

% pressure
p0 = 4.5e-3; % [GPa]
BC.initP = p0 * ones(MeshP.nDOF,1);

%% Dirichlet BCs - solid
% displacement u=0 at boundaries
BC.fixed_u = [MeshU.bottom_dofy; MeshU.top_dofy; MeshU.left_dofx; MeshU.right_dofx];
% fixed DOF values
BC.fixed_u_value = zeros(length(BC.fixed_u),1);
% free nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid
% nodes at injection well 1
MeshP.nodesInj1 = [5; 71; 72; 73; 74; 75; 76; 77];
% nodes at injection well 2
MeshP.nodesInj2 = [6; 78; 79; 80; 81; 82; 83; 84];
% pressure at injection wells [GPa]
pInj = 3e-3; 
% fixed DOFs (nodes at injection wells)
BC.fixed_p = [MeshP.nodesInj1; MeshP.nodesInj2];
% fixed DOF values
BC.fixed_p_value = pInj * ones(length(BC.fixed_p),1);
% free nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% prescribed traction [GN/m2]
BC.tractionNodes = [];

% prescribed traction
BC.traction = 10e-3;
% traction interpolation (needed for traction applied in wells); 1 - true, 0 - false
BC.tractionInterp = 1;
% traction unit normal in local coords
normal_L = [0; -1];

% nodes at wells
BC.tractionNodes1 = [5; 139; 140; 141; 142; 143; 144; 145; 146; 147; 148; 149; 150; 151; 152; 153];
BC.tractionNodes2 = [6; 154; 155; 156; 157; 158; 159; 160; 161; 162; 163; 164; 165; 166; 167; 168];
BC.tractionNodes = [BC.tractionNodes1; BC.tractionNodes2];

% traction force vector
BC.tractionForce1 = zeros(length(BC.tractionNodes1),2);
BC.tractionForce2 = zeros(length(BC.tractionNodes2),2);

% coordinates of well centre
centre1 = [60, 15];
centre2 = [60, 10];

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
    BC.tractionForce1(i,:) = BC.traction * normal_G';
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
    BC.tractionForce2(i,:) = BC.traction * normal_G';
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
% BC.b = @(x)[-10e-3; -10e-3];  
BC.b = @(x)[];

%% Neumann BCs - fluid
% distributed flux [m3/s]
% impervious at bottom, left, and right
BC.fluxNodes = [];
BC.fluxValue = zeros(length(BC.fluxNodes),1);

% % flux at wells
% BC.flux = 1;
% 
% % well 1
% BC.fluxNodes1 = [5; 71; 72; 73; 74; 75; 76; 77];
% % find surface area
% surfs1 = zeros(length(BC.fluxNodes1), 2);
% L = 0;
% count = 1;
% for e = 1:MeshP.ne
%     conn_e = MeshP.conn(e,:);
%     for m = 1:length(BC.fluxNodes1)
%         nodeM = BC.fluxNodes1(m);
%         for n = 1:length(BC.fluxNodes1)
%             nodeN = BC.fluxNodes1(n);
%             if m ~= n && any(conn_e == nodeM) && any(conn_e == nodeN) && (isequal([nodeM;nodeN], [conn_e(1);conn_e(2)]) || isequal([nodeM;nodeN], [conn_e(2);conn_e(3)]) || isequal([nodeM;nodeN], [conn_e(3);conn_e(4)]) || isequal([nodeM;nodeN], [conn_e(4);conn_e(1)]))
%                 surfs1(count,1) = nodeM;
%                 surfs1(count,2) = nodeN;
%                 L = L + sqrt((MeshP.coords(nodeM,1) - MeshP.coords(nodeN,1))^2 +...
%                     (MeshP.coords(nodeM,2) - MeshP.coords(nodeN,2))^2);
%                 count = count + 1;
%             end
%         end
%     end
% end
% % total flux
% Flux1 = BC.flux * L/length(BC.fluxNodes1);
% BC.fluxValue1 = Flux1*ones(size(BC.fluxNodes1));
% 
% % well 2
% BC.fluxNodes2 = [6; 78; 79; 80; 81; 82; 83; 84];
% % find surface area
% surfs2 = zeros(length(BC.fluxNodes2), 2);
% L2 = 0;
% count = 1;
% for e = 1:MeshP.ne
%     conn_e = MeshP.conn(e,:);
%     for m = 1:length(BC.fluxNodes2)
%         nodeM = BC.fluxNodes2(m);
%         for n = 1:length(BC.fluxNodes2)
%             nodeN = BC.fluxNodes2(n);
%             if m ~= n && any(conn_e == nodeM) && any(conn_e == nodeN) && (isequal([nodeM;nodeN], [conn_e(1);conn_e(2)]) || isequal([nodeM;nodeN], [conn_e(2);conn_e(3)]) || isequal([nodeM;nodeN], [conn_e(3);conn_e(4)]) || isequal([nodeM;nodeN], [conn_e(4);conn_e(1)]))
%                 surfs2(count,1) = nodeM;
%                 surfs2(count,2) = nodeN;
%                 L2 = L2 + sqrt((MeshP.coords(nodeM,1) - MeshP.coords(nodeN,1))^2 +...
%                     (MeshP.coords(nodeM,2) - MeshP.coords(nodeN,2))^2);
%                 count = count + 1;
%             end
%         end
%     end
% end
% % total flux
% Flux2 = BC.flux * L2/length(BC.fluxNodes2);
% BC.fluxValue2 = Flux2*ones(size(BC.fluxNodes2));
% 
% % both wells
% BC.fluxNodes = [BC.fluxNodes1; BC.fluxNodes2];
% BC.fluxValue = [BC.fluxValue1; BC.fluxValue2];

% point flux [m/s]
BC.pointFlux = [];

% flux source [m3/s/m3]
BC.s = @(x)[]; 

%% Porosity BCs
if ~Control.Biotmodel
    BC.fixed_n = [];
    BC.free_n = setdiff(MeshN.DOF, BC.fixed_n);
    BC.fixed_n_value = zeros(length(BC.fixed_n),1);
end

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
Control.dt = 1;  % time step [s]
Control.tend = 100;   % final simulation time [s]

Control.beta = 1; % beta-method time discretization -- beta = 1 Backward Euler; beta = 0.5 Crank-Nicolson

Control.plotu = 223*2; 
Control.plotp = 139;

% plot analytical solution (valid for 1D problems with Material.Minv == 0)
Control.plotansol = 0; % 1 = true; 0 = false

% solve in the frequency domain
Control.freqDomain = 0;  % 1 = true; 0 = false

end