% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function ComputeTimeStepError_UP()
% Compute error using the L2 norm of the time discretization

clearvars
clear, clc

%% Initialize variables
nsims = 4; % number of simulations
h = zeros(nsims,1); % time step

% discretization of time domain
nsd = 1; % number of dimensions
type = 'L2'; % element type
field = 'p'; % equivalent to scalar field

% The name of the files need to be updated according to the desired plots
%% Mesh 1
% load file
load Results_m76ManUPUSpanosdt.mat
L = Control.tend; % final simulation time
% discretization of time domain
ne = round(Control.tend/Control.dt); % number of elements
% generate mesh
Mesh1 = Build1DMesh(nsd, ne, L, type, field);
% pressure values
dp1 = Plot.p_time;
% solid displacement values
dus1 = Plot.u_time;
% time step
h(1) = Control.dt;

%% Mesh 2
load Results_m77ManUPUSpanosdt.mat
% discretization of time domain
ne = round(Control.tend/Control.dt); % number of elements
% generate mesh
Mesh2 = Build1DMesh(nsd, ne, L, type, field);
% pressure values
dp2 = Plot.p_time;
% solid displacement values
dus2 = Plot.u_time;
% time step
h(2) = Control.dt;

%% Mesh 3
load Results_m78ManUPUSpanosdt.mat
% discretization of time domain
ne = round(Control.tend/Control.dt); % number of elements
% generate mesh
Mesh3 = Build1DMesh(nsd, ne, L, type, field);
% pressure values
dp3 = Plot.p_time;
% solid displacement values
dus3 = Plot.u_time;
% time step
h(3) = Control.dt;

%% Mesh 4
load Results_m79ManUPUSpanosdt.mat
% discretization of time domain
ne = round(Control.tend/Control.dt); % number of elements
% generate mesh
Mesh4 = Build1DMesh(nsd, ne, L, type, field);
% pressure values
dp4 = Plot.p_time;
% solid displacement values
dus4 = Plot.u_time;
% time step
h(4) = Control.dt;

%% Reference mesh in time
% load ref results
load Results_m80ManUPUSpanosdt.mat
% discretization of time domain
nsd = 1; % number of dimensions
% discretization of time domain
ne = round(Control.tend/Control.dt); % number of elements
% generate mesh
Mesh_Exact = Build1DMesh(nsd, ne, L, type, field);
% pressure values
dp_exact = Plot.p_time;
% solid displacement values
dus_exact = Plot.u_time;

%% Beta parameter
beta = Control.beta;

%% Compute L2 error norm
% initialize variables
eL2p = zeros(nsims,1);
eL2us = zeros(nsims,1);

% loop over meshes
for sim = 1:nsims
    switch sim
        case 1
            dp = dp1;
            dus = dus1;
            Mesh = Mesh1;
        case 2
            dp = dp2;
            dus = dus2;
            Mesh = Mesh2;
        case 3
            dp = dp3;
            dus = dus3;
            Mesh = Mesh3;
        case 4
            dp = dp4;
            dus = dus4;
            Mesh = Mesh4;
        case 5
            dp = d5;
            dus = dus5;
            Mesh = Mesh5;
        case 6
            dp = d6;
            dus = dus6;
            Mesh = Mesh6;
    end
        
    % Calculate error norm
    eL2p_num = 0;
    eL2p_den = 0;
    
    eL2us_num = 0;
    eL2us_den = 0;
    
    % higher quadrature definition
    Control.nqP = 8;
    Quad = GlobalQuad(Mesh, Control);

    % interpolate 'approx.' mesh at 'exact' mesh nodes
    dp_interp = zeros(Mesh_Exact.nn,1);
    dus_interp = zeros(Mesh_Exact.nn,1);
    
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
               dp_ap = dp(conn_ap);
               % nodal values (solid displacement)
               dus_ap = dus(conn_ap);
               % exact mesh node in parent coords
               coord_ex_parent = (2*coord_ex-coords_ap(2)-coords_ap(1))/(coords_ap(2)-coords_ap(1));
               % approx mesh element shape functions
               N = lagrange_basis(Mesh, coord_ex_parent);
               dp_interp(i) = N'*dp_ap;
               dus_interp(i) = N'*dus_ap;
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
            uph = N*dp_interp(Mesh_Exact.xdofs(conn_e));
            % exact displacement at quadrature point
            upe = N*dp_exact(Mesh_Exact.xdofs(conn_e));
            
            % approximated displacement at quadrature point
            ush = N*dus_interp(Mesh_Exact.xdofs(conn_e));
            % exact displacement at quadrature point
            use = N*dus_exact(Mesh_Exact.xdofs(conn_e));
            
            % L2 norm
            eL2p_num = eL2p_num + ((uph - upe).^2) * Quad.w(ip,1) * Jdet;
            eL2p_den = eL2p_den + (upe.^2) * Quad.w(ip,1) * Jdet;
            
            eL2us_num = eL2us_num + ((ush - use).^2) * Quad.w(ip,1) * Jdet;
            eL2us_den = eL2us_den + (use.^2) * Quad.w(ip,1) * Jdet;
        end
    end
    
    eL2p(sim) = sqrt(eL2p_num/eL2p_den);
    eL2us(sim) = sqrt(eL2us_num/eL2us_den);
end

% Determine slope of L2 norm
pL2p = polyfit(log(h), log(eL2p),1);
m_L2p = pL2p(1);

pL2us = polyfit(log(h), log(eL2us),1);
m_L2us = pL2us(1);

figure;
loglog(h,eL2p,'g*', 'LineWidth', 1.5);
hold on
loglog(h,exp(pL2p(2))*h.^pL2p(1),'g', 'LineWidth', 1.5);
grid on
xlabel('Mesh size (m)');
ylabel('L2-norm');
title(sprintf('Time - p Convergence beta = %.2f, order %.4f', beta, m_L2p));
hold off

figure;
loglog(h,eL2us,'m*', 'LineWidth', 1.5);
hold on
loglog(h,exp(pL2us(2))*h.^pL2us(1),'m', 'LineWidth', 1.5);
grid on
xlabel('Mesh size (m)');
ylabel('L2-norm');
title(sprintf('Time - us Convergence beta = %.2f, order %.4f', beta, m_L2us));
hold off

end