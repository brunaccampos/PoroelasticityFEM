function [Material, MeshU, MeshP, MeshN, BC, Control] = PlateWithHole_Elasticity(config_dir, progress_on)
% ------------------------------------------------------------------------
% 2D elasticity problem - steady case
% Adapted from: https://github.com/GCMLab
% ------------------------------------------------------------------------

%% Material properties
% elasticity modulus [Pa]
Material.E = 2e11;
% Poisson's ratio
Material.nu = 0.3;

% porous media permeability [m2/Pa s]
Material.kf = 0;
% Biot's coefficient
Material.alpha = 0;
% 1/Q (related to storage coefficient)
Material.Minv = 0;
% poroelasticity model
Control.Biotmodel = 1;

% lumped mass matrix - 0: false, 1: true
Material.lumpedMass = 1;

% thickness 
% 1D: cross sectional area [m2]
% 2D: out of plane thickness [m]
Material.t = 1;

% constititive law - 'PlaneStress' or 'PlaneStrain'
% Note: use 'PlaneStrain' for 1D or 2D poroelasticity
Material.constLaw = 'PlaneStrain';

%% Mesh Properties
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
        meshFileNameU = 'Mesh Files\PlateWithHoleQ4.msh';
        MeshU = BuildMesh_GMSH(meshFileNameU, fieldU, nsd, config_dir, progress_on);
        %%%% pressure field
        fieldP = 'p';
        meshFileNameP = 'Mesh Files\PlateWithHoleQ4.msh';
        MeshP = BuildMesh_GMSH(meshFileNameP, fieldP, nsd, config_dir, progress_on);
        %%%% porosity field
        if ~Control.Biotmodel
            fieldN = 'n';
            meshFileNameN = 'Mesh Files\PlateWithHoleQ4.msh';
            MeshN = BuildMesh_GMSH(meshFileNameN, fieldN, nsd, config_dir, progress_on);
        else
            MeshN = [];
        end
end

%% Dirichlet BCs - solid
% column vector of prescribed displacement dof
BC.fixed_u = [MeshU.left_dofx; MeshU.bottom_dofy];
% prescribed displacement for each dof [u1; u2; ...] [m]
BC.fixed_u_value = zeros(length(BC.fixed_u),1);
% free nodes
BC.free_u = setdiff(MeshU.DOF, BC.fixed_u);

%% Dirichlet BCs - fluid
BC.fixed_p = 1:MeshP.nDOF;
% fixed DOF values
BC.fixed_p_value = zeros(length(BC.fixed_p),1);
% free nodes
BC.free_p = setdiff(MeshP.DOF, BC.fixed_p);

%% Neumann BCs - solid
% column vector of prescribed traction nodes
BC.tractionNodes = MeshU.right_nodes;
% distributed traction [N/m2]
BC.traction = 10e3; % uniform tensile stress applied to right edge
Fright = BC.traction*max(MeshU.coords(:,2))/((length(MeshU.right_nodes) - 1)/2);
BC.tractionForce = zeros(length(BC.tractionNodes),2);
switch MeshU.type
    case 'Q4'
        BC.tractionForce = [Fright/2 *ones(size(MeshU.right_nodes)),     zeros(size(MeshU.right_nodes))];
    case 'Q9'
        for n = 1:length(BC.tractionNodes)
            if any( BC.tractionNodes(n) == MeshU.conn(:,1:4),'all') % then node is a corner node
                BC.tractionForce(n,:) = [Fright/3, 0];
            else % then node is a midside node
                BC.tractionForce(n,:) = [Fright*2/3,0];
            end
        end
end

% find the nodes in the top left and bottom right corners
botrightnode = find(MeshU.coords(BC.tractionNodes,2) == min(MeshU.coords(:,2)));
toprightnode  = find(MeshU.coords(BC.tractionNodes,2) == max(MeshU.coords(:,2)));

BC.tractionForce(botrightnode,1) = BC.tractionForce(botrightnode,1)/2;
BC.tractionForce(toprightnode,1) = BC.tractionForce(toprightnode,1)/2;

% point loads [N]
BC.pointLoad = [];

% body force [N/m3]
BC.b = @(x)[];  

%% Neumann BCs - fluid
% distributed flux [m/s]
BC.fluxNodes = 1:MeshP.nDOF;
BC.fluxValue = zeros(length(BC.fluxNodes),1);

% point flux [m3/s]
BC.pointFlux = [];

% flux source [m3/s/m3]
BC.s = @(x)[]; 

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
Control.dt = 1;  % time step
Control.tend = 1;   % final simulation time

Control.beta = 1; % beta-method time discretization -- beta = 1 Backward Euler; beta = 0.5 Crank-Nicolson

Control.plotu = 2;
Control.plotp = 2;

% plot analytical solution (valid for 1D problems with Material.Minv == 0)
Control.plotansol = 0; % 1 = true; 0 = false

% solve in the frequency domain
Control.freqDomain = 0;  % 1 = true; 0 = false

end