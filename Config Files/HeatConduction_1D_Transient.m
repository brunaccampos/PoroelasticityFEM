function [Material, MeshU, MeshP, MeshN, BC, Control] = HeatConduction_1D_Transient(config_dir, progress_on)
% 
% Configuration File
% ------------------------------------------------------------------------

%% Material properties
% porous media permeability [m2/Pa s]
Material.kf = 0.002;
% 1/Q (related to storage coefficient)
Material.Minv = 1;
% Poroelasticity model
Control.PMmodel = 'Tr1_Biot_UP';

% elasticity modulus [Pa]
Material.E = 0;
% Poisson's ratio
Material.nu = 0;
% Biot's coefficient
Material.alpha = 0;

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
MeshType = 'Manual';

switch MeshType
    case 'Manual'
        % number of space dimensions
        nsd = 1;
        % number of elements
        ne = 8;
        % column size [m]
        L = 2;
        %%%% solid displacement field
        typeU = 'L2';
        fieldU = 'u';
        MeshU = Build1DMesh(nsd, ne, L, typeU, fieldU);
        %%%% fluid pressure field
        typeP = 'L3';
        fieldP = 'p';
        MeshP = Build1DMesh(nsd, ne, L, typeP, fieldP);
        %%%% porosity field
        if contains(Control.PMmodel, 'UPN')
            typeN = 'L2';
            fieldN = 'n';
            MeshN = Build1DMesh(nsd, ne, L, typeN, fieldN);
        else
            MeshN = [];
        end
    case 'Gmsh'
        % Version 2 ASCII
        % number of space dimensions
        nsd = 2;
        %%%% displacement field
        fieldU = 'u';
        meshFileNameU = 'Column2DQ9.msh';
        MeshU = BuildMesh_GMSH(meshFileNameU, fieldU, nsd, config_dir, progress_on);
        %%%% pressure field
        fieldP = 'p';
        meshFileNameP = 'Column2DQ4.msh';
        MeshP = BuildMesh_GMSH(meshFileNameP, fieldP, nsd, config_dir, progress_on);
        %%%% porosity field
        if contains(Control.PMmodel, 'UPN')
            fieldN = 'n';
            meshFileNameN = 'Column2DQ4.msh';
            MeshN = BuildMesh_GMSH(meshFileNameN, fieldN, nsd, config_dir, progress_on);
        else
            MeshN = [];
        end
end

%% Initial conditions
% displacement
BC.initU = [];

% pressure
BC.initP = 50*ones(MeshP.nDOF,1);

%% Find nodes for prescribed BCs
% find top and bottom nodes for displacement field
BC.top_node_u = find(MeshU.coords == max(MeshU.coords));
BC.bottom_node_u = find(MeshU.coords == min(MeshU.coords));

% find top and bottom nodes for pressure field
BC.top_node_p = find(MeshP.coords == max(MeshP.coords));
BC.bottom_node_p = find(MeshP.coords == min(MeshP.coords));

%% Dirichlet BCs - solid
% column vector of prescribed displacement dof
BC.fixed_u = 1:MeshU.nDOF;
% prescribed displacement for each dof [u1; u2; ...] [m]
BC.fixed_u_value = zeros(length(BC.fixed_u),1);
% free nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid
%   pressure p=0 at the top
BC.fixed_p = [BC.bottom_node_p; BC.top_node_p];
% fixed DOF values
BC.fixed_p_value = [100; 0];
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
% point flux [m3/s]
BC.pointFlux = [];

% distributed flux [m/s]
BC.fluxNodes = [];

% flux source [m3/s/m3]
BC.s = @(x)[]; 

%% Porosity BCs
if contains(Control.PMmodel, 'UPN')
    BC.fixed_n = [];
    BC.free_n = setdiff(MeshN.DOF, BC.fixed_n);
    BC.fixed_n_value = zeros(length(BC.fixed_n),1);
end

%% Quadrature order
Control.nqU = 2;
Control.nqP = 2;

%% Solution parameters
% tag used for computing analytical solution
% 1 = uncoupled problem (elasticity, heat transfer, etc)
% 0 = coupled problem (Biot, Spanos model)
Control.uncoupled = 1; 

% basic time step controls
Control.dt = 1e-2;  % time step
Control.tend = 10;   % final simulation time

Control.beta = 1; % beta-method time discretization -- beta = 1 Backward Euler; beta = 0.5 Crank-Nicolson

% DOF to plot graphs
Control.plotu = 1;
Control.plotp = 1;

% plot analytical solution (valid for 1D problems with Material.Minv == 0)
Control.plotansol = 1; % 1 = true; 0 = false

% solve in the frequency domain
Control.freqDomain = 0;  % 1 = true; 0 = false

%% Analytical solution
Control.uan_symb = @(x) x*0;
Control.u_an = Control.uan_symb(MeshU.coords);

aux=0;
syms x
N=1000;
for k=1:N
    aux = aux + (1/k)*exp(-Material.kf*k^2*pi()^2*Control.tend)*sin(k*pi()*x);
end
steady = (BC.fixed_p_value(2) - BC.fixed_p_value(1))*x/max(MeshP.coords) + BC.fixed_p_value(1);

Control.pan_symb = steady - (BC.fixed_p_value(1)/pi()) * aux; % transient solution


aux=0;
x = MeshP.coords;
N=1000;
for k=1:N
    aux = aux + (1/k)*exp(-Material.kf*k^2*pi()^2*Control.tend)*sin(k*pi()*x);
end
steady = (BC.fixed_p_value(2) - BC.fixed_p_value(1))*x/max(MeshP.coords) + BC.fixed_p_value(1);

Control.p_an = steady - (BC.fixed_p_value(1)/pi()) * aux; % transient solution

end
