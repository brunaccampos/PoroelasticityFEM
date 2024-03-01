function [Material, MeshU, MeshP, MeshN, BC, Control] = PatchTestD(config_dir, progress_on, ~, ~)
% ------------------------------------------------------------------------
% All nodes are restrained and nodal pressure values
% are specfied according to the exact solution. The error between the FEA
% and exact solutions is then calculated. The FEA approximate solution
% should be exact.
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

%% Material properties
% diffusivity coefficient [m2/s]
Material.kf = 0.139e-4;
% elasticity modulus [Pa]
Material.E = 0;
% Poisson's ratio
Material.nu = 0;
% 1/Q (related to storage coefficient)
Material.Minv = 0;
% Biot's coefficient
Material.alpha = 0;

% thickness 
% 1D: cross sectional area [m2]
% 2D: out of plane thickness [m]
Material.t = 1;

% constititive law - 'PlaneStress' or 'PlaneStrain'
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
        L = 1;
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
        meshFileNameU = 'Mesh Files\PatchTest.msh';
        MeshU = BuildMesh_GMSH(meshFileNameU, fieldU, nsd, config_dir, progress_on);
        %%%% pressure field
        fieldP = 'p';
        meshFileNameP = 'Mesh Files\PatchTest.msh';
        MeshP = BuildMesh_GMSH(meshFileNameP, fieldP, nsd, config_dir, progress_on);
        %%%% porosity field
        if contains(Control.PMmodel, 'UPN')
            fieldN = 'n';
            meshFileNameN = 'Mesh Files\PatchTest.msh';
            MeshN = BuildMesh_GMSH(meshFileNameN, fieldN, nsd, config_dir, progress_on);
        else
            MeshN = [];
        end
end

% flux [m/s] (same in both directions)
BC.flux = 1;

%% Dirichlet BCs - solid
% prescribed displacement
BC.fixed_u = 1:MeshU.nDOF;
BC.fixed_u_value = @(t) zeros(length(BC.fixed_u),1);
% free displacement nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid
% pressure according to exact solution
BC.p = @(x) -1/Material.kf * BC.flux * (x(:,1) + x(:,2));
% fixed nodes
BC.fixed_p = 1:MeshP.nDOF;
% prescribed pressure
BC.fixed_p_value = @(t) BC.p(MeshP.coords(BC.fixed_p, :));
% free pressure nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% point load [N]
BC.pointLoadValue = 0;
BC.pointLoadNodes = 1:MeshU.nDOF;
BC.pointLoad = zeros(MeshU.nDOF,1);
BC.pointLoad(BC.pointLoadNodes) = BC.pointLoadValue;

% traction interpolation (needed for traction applied in wells); 1 - true, 0 - false
BC.tractionInterp = 0;

% distributed load [N/m]
BC.tractionNodes = [];

% body force
BC.b = @(x,t)[];

%% Neumann BCs - fluid
BC.fluxNodes = [];

% point flux [m/s]
BC.pointFlux = [];

% flux source
BC.s = @(x,t)[]; 

%% Quadrature order
Control.nqU = 2;
Control.nqP = 2;

%% Time step controls
Control.dt = 1;  % time step
Control.tend = 1; % final simulation time

% Beta method
% beta = 1 Backward Euler; beta = 0.5 Crank-Nicolson
Control.beta = 1; 

%% Plot data
% DOF to plot graphs
Control.plotu = 1;
Control.plotp = 1;

end