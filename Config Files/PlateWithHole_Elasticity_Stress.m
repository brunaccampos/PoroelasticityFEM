function [Material, MeshU, MeshP, MeshN, BC, Control] = PlateWithHole_Elasticity_Stress(config_dir, progress_on)
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
Control.PMmodel = 'Tr1_Biot_UP';

% lumped mass matrix - 0: false, 1: true
Material.lumpedMass = 1;

% thickness 
% 1D: cross sectional area [m2]
% 2D: out of plane thickness [m]
Material.t = 1;

% constititive law - 'PlaneStress' or 'PlaneStrain'
% Note: use 'PlaneStrain' for 1D or 2D poroelasticity
Material.constLaw = 'PlaneStress';

%% In situ stress field
% [GPa] negative for compression
BC.S0 = [];

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
        meshFileNameU = 'Mesh Files\PlateWithHoleQ4.msh';
        MeshU = BuildMesh_GMSH(meshFileNameU, fieldU, nsd, config_dir, progress_on);
        %%%% pressure field
        fieldP = 'p';
        meshFileNameP = 'Mesh Files\PlateWithHoleQ4.msh';
        MeshP = BuildMesh_GMSH(meshFileNameP, fieldP, nsd, config_dir, progress_on);
        %%%% porosity field
        if contains(Control.PMmodel, 'UPN')
            fieldN = 'n';
            meshFileNameN = 'Mesh Files\PlateWithHoleQ4.msh';
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
% column vector of prescribed displacement dof
BC.fixed_u = [];
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
BC.tractionNodes = [1; 5; 6; 7; 8; 9; 10; 11; 12];

% distributed traction [N/m2]
BC.traction = 10e-3;
% traction unit normal in local coords
normal_L = [0; -1];
% traction force vector
BC.tractionForce = zeros(length(BC.tractionNodes),2);

% find surfaces with traction applied
surfsTrac = zeros(length(BC.tractionNodes), 2);
count = 1;
for e = 1:MeshU.ne
    conn_e = MeshU.conn(e,:);
    for m = 1:length(BC.tractionNodes)
        nodeM = BC.tractionNodes(m);
        for n = 1:length(BC.tractionNodes)
            nodeN = BC.tractionNodes(n);
            if m ~= n && any(conn_e == nodeM) && any(conn_e == nodeN) && (isequal([nodeM;nodeN], [conn_e(1);conn_e(2)]) || isequal([nodeM;nodeN], [conn_e(2);conn_e(3)]) || isequal([nodeM;nodeN], [conn_e(3);conn_e(4)]) || isequal([nodeM;nodeN], [conn_e(4);conn_e(1)]))
                surfsTrac(count,:) = [nodeM, nodeN];
                % surface nodes
                node1 = surfsTrac(count,1);
                node2 = surfsTrac(count,2);
                % length of each segment
                L = sqrt((MeshU.coords(node1,1) - MeshU.coords(node2,1))^2 +...
                    (MeshU.coords(node1,2) - MeshU.coords(node2,2))^2);
                % sin and cos for segment
                cos12 = abs(MeshU.coords(node1,2) - MeshU.coords(node2,2)) / L;
                sin12 = abs(MeshU.coords(node1,1) - MeshU.coords(node2,1)) / L;
                % rotation matrices for normal vector
                rot12 = [cos12, -sin12; sin12, cos12];
                % equivalent forces
                force12 = rot12 * normal_L * BC.traction * L / ((length(BC.tractionNodes)-1)/2);
                % store in nodes
                index1 = find(BC.tractionNodes == node1);
                index2 = find(BC.tractionNodes == node2);
                BC.tractionForce(index1,:) = BC.tractionForce(index1,:) + force12';
                BC.tractionForce(index2,:) = BC.tractionForce(index2,:) + force12';

                % update counter
                count = count + 1;
            end
        end
    end
end

BC.tractionForce = BC.tractionForce/2;
indexBC1 = find(BC.tractionNodes == 1);
indexBC2 = find(BC.tractionNodes == 5);

BC.tractionForce(indexBC1,:) = BC.tractionForce(indexBC1,:)/2;
BC.tractionForce(indexBC2,:) = BC.tractionForce(indexBC2,:)/2;


% Fright = BC.traction*max(MeshU.coords(:,2))/((length(MeshU.right_nodes) - 1)/2);
% BC.tractionForce = zeros(length(BC.tractionNodes),2);
% switch MeshU.type
%     case 'Q4'
%         BC.tractionForce = [Fright/2 *ones(size(MeshU.right_nodes)),     zeros(size(MeshU.right_nodes))];
%     case 'Q9'
%         for n = 1:length(BC.tractionNodes)
%             if any( BC.tractionNodes(n) == MeshU.conn(:,1:4),'all') % then node is a corner node
%                 BC.tractionForce(n,:) = [Fright/3, 0];
%             else % then node is a midside node
%                 BC.tractionForce(n,:) = [Fright*2/3,0];
%             end
%         end
% end
% 
% % find the nodes in the top left and bottom right corners
% botrightnode = find(MeshU.coords(BC.tractionNodes,2) == min(MeshU.coords(:,2)));
% toprightnode  = find(MeshU.coords(BC.tractionNodes,2) == max(MeshU.coords(:,2)));
% 
% BC.tractionForce(botrightnode,1) = BC.tractionForce(botrightnode,1)/2;
% BC.tractionForce(toprightnode,1) = BC.tractionForce(toprightnode,1)/2;

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

%% Solution parameters
% tag used for computing analytical solution
% 1 = uncoupled problem (elasticity, heat transfer, etc)
% 0 = coupled problem (Biot, Spanos model)
Control.uncoupled = 0; 

% basic time step controls
Control.dt = 1;  % time step
Control.tend = 1;   % final simulation time

Control.beta = 1; % beta-method time discretization -- beta = 1 Backward Euler; beta = 0.5 Crank-Nicolson

% DOF to plot graphs
Control.plotu = 2;
Control.plotp = 2;

% plot analytical solution (valid for 1D problems with Material.Minv == 0)
Control.plotansol = 0; % 1 = true; 0 = false

% solve in the frequency domain
Control.freqDomain = 0;  % 1 = true; 0 = false

end