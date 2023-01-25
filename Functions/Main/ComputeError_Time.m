function ComputeError_Time()
% ------------------------------------------------------------------------
% Compute error using the Le norm of the displacements
% ------------------------------------------------------------------------

clearvars
clear, clc
close all

%% Initialize variables
nsims = 4; % number of simulations
h = zeros(nsims,1); % time step

% discretization of time domain
nsd = 1; % number of dimensions
type = 'L2'; % element type
field = 'p'; % equivalent to scalar field

%% Mesh 1
% load file
load Results_m1_time_Boone1e-1_beta05.mat

L = Control.tend; % final simulation time
beta = Control.beta; % beta method parameter

% discretization of time domain
ne = round(Control.tend/Control.dt); % number of elements
% generate mesh
Mesh1 = Build1DMesh(nsd, ne, L, type, field);
% pressure values
d1 = Plot.p_time;
% time step
h(1) = Control.dt;

%% Mesh 2
load Results_m2_time_Boone1e-2_beta05.mat

% discretization of time domain
ne = round(Control.tend/Control.dt); % number of elements
% generate mesh
Mesh2 = Build1DMesh(nsd, ne, L, type, field);
% pressure values
d2 = Plot.p_time;
% time step
h(2) = Control.dt;

%% Mesh 3
load Results_m3_time_Boone1e-3_beta05.mat

% discretization of time domain
ne = round(Control.tend/Control.dt); % number of elements
% generate mesh
Mesh3 = Build1DMesh(nsd, ne, L, type, field);
% pressure values
d3 = Plot.p_time;
% time step
h(3) = Control.dt;

%% Mesh 4
load Results_m4_time_Boone1e-4_beta05.mat

% discretization of time domain
ne = round(Control.tend/Control.dt); % number of elements
% generate mesh
Mesh4 = Build1DMesh(nsd, ne, L, type, field);
% pressure values
d4 = Plot.p_time;
% time step
h(4) = Control.dt;

%% Reference mesh in time
% load ref results
load Results_m5_time_Boone1e-5_beta05.mat

% discretization of time domain
nsd = 1; % number of dimensions
L = Control.tend; % final simulation time
ne = round(Control.tend/Control.dt); % number of elements
type = 'L2'; % element type
field = 'p'; %scalar field

% pressure values
d_exact = Plot.p_time;

% generate mesh
Mesh_Exact = Build1DMesh(nsd, ne, L, type, field);

%% Compute L2 error norm
% initialize variables
eL2 = zeros(nsims,1);

% loop over meshes
for sim = 1:nsims
    switch sim
        case 1
            d = d1;
            Mesh = Mesh1;
        case 2
            d = d2;
            Mesh = Mesh2;
        case 3
            d = d3;
            Mesh = Mesh3;
        case 4
            d = d4;
            Mesh = Mesh4;
        case 5
            d = d5;
            Mesh = Mesh5;
        case 6
            d = d6;
            Mesh = Mesh6;
    end
        
    % Calculate error norm
    eL2_num = 0;
    eL2_den = 0;
    
    % higher quadrature definition
    Control.nqU = 8;
    MeshU.field = 'u';
    Quad = GlobalQuad(MeshU, Control);

    % interpolate 'approx.' mesh at 'exact' mesh nodes
    d_interp = zeros(Mesh_Exact.nn,1);
    
    % loop over exact mesh nodes
    for i = 1:Mesh_Exact.nn
       coord_ex = Mesh_Exact.coords(i);
       % loop over approx mesh elements
       for e = 1:Mesh.ne
           % element connectivity
           conn_ap = Mesh.conn(e,:);
           % element nodes coordinates
           coords_ap = Mesh.coords(conn_ap);
           % check in which element of the approx mesh the exact mesh node is located
           if coords_ap(1) <= coord_ex && coord_ex <= coords_ap(2)
               % nodal values (pressure)
               d_ap = d(conn_ap);
               % exact mesh node in parent coords
               coord_ex_parent = (2*coord_ex-coords_ap(2)-coords_ap(1))/(coords_ap(2)-coords_ap(1));
               % approx mesh element shape functions
               N = lagrange_basis(Mesh, coord_ex_parent);
               d_interp(i) = N'*d_ap;
           end
       end
    end
    
    % loop over elements
    for i = 1:Mesh_Exact.ne
        % element connectivity
        conn_e = Mesh_Exact.conn(i,:);
        % global coordinates
        gcoords = Mesh_Exact.coords(conn_e,:);
                
        % loop over integration points
        for ip = 1:Quad.nq
            % N matrices
            N = getN(Mesh_Exact, Quad, ip);
            % N derivatives
            dN = getdN(Mesh_Exact, Quad, ip);
            % Jacobian matrix
            J = dN*gcoords;
            % Jacobian determinant
            Jdet = det(J);
            
            % approximated displacement at quadrature point
            uh = N*d_interp(Mesh_Exact.xdofs(conn_e));
            % exact displacement at quadrature point
            ue = N*d_exact(Mesh_Exact.xdofs(conn_e));
            
            % L2 norm
            eL2_num = eL2_num + ((uh - ue).^2) * Quad.w(ip,1) * Jdet;
            eL2_den = eL2_den + (ue.^2) * Quad.w(ip,1) * Jdet;
        end
    end
    
    eL2(sim) = sqrt(eL2_num/eL2_den);
end

% Determine slope of L2 norm
pL2 = polyfit(log(h), log(eL2),1);
m_L2 = pL2(1);

figure;
loglog(h,eL2,'g-o', 'LineWidth', 1.5);
xlabel('Mesh size (m)');
ylabel('L2-norm');
title(sprintf('Convergence for time beta = %.2f, order %.4f', beta, m_L2));



end