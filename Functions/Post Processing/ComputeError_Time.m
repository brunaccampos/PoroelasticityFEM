function ComputeError_Time()
% ------------------------------------------------------------------------
% Compute error using the Le norm of the displacements
% ------------------------------------------------------------------------

clearvars
clear, clc
close all

%% Number of sims
nsims = 5;

nsd = 1;
L = 50; % t end
type = 'L2';

%% Mesh 1
% load file
load Results_m1Boone0.1Spanos_Beta1.mat
% mesh
ne1 = L/Control.dt; % dt
Mesh1 = Build1DMesh(nsd, ne1, L, type);
% FEM approximation
d1 = Plot.u_time;
% exact solution
d_exact1 = Plot.uan_time;

%% Mesh 2
load Results_m2Boone0.05Spanos_Beta1.mat
% mesh
ne2 = L/Control.dt; % dt
Mesh2 = Build1DMesh(nsd, ne2, L, type);
% FEM approximation
d2 = Plot.u_time;
% exact solution
d_exact2 = Plot.uan_time;

%% Mesh 3
load Results_m3Boone0.025Spanos_Beta1.mat
% mesh
ne3 = L/Control.dt; % dt
Mesh3 = Build1DMesh(nsd, ne3, L, type);
% FEM approximation
d3 = Plot.u_time;
% exact solution
d_exact3 = Plot.uan_time;

%% Mesh 4
load Results_m4Boone0.0125Spanos_Beta1.mat
% mesh
ne4 = L/Control.dt; % dt
Mesh4 = Build1DMesh(nsd, ne4, L, type);
% FEM approximation
d4 = Plot.u_time;
% exact solution
d_exact4 = Plot.uan_time;

%% Mesh 5
load Results_m5Boone0.00625Spanos_Beta1.mat
% mesh
ne5 = L/Control.dt; % dt
Mesh5 = Build1DMesh(nsd, ne5, L, type);
% FEM approximation
d5 = Plot.u_time;
% exact solution
d_exact5 = Plot.uan_time;

%% Compute L2 error norm
% initialize variables
eL2 = zeros(nsims,1);
h = zeros(nsims,1);

% loop over meshes
for sim = 1:nsims
    switch sim
        case 1
            d = d1;
            d_exact = d_exact1;
            Mesh = Mesh1;
        case 2
            d = d2;
            d_exact = d_exact2;
            Mesh = Mesh2;
        case 3
            d = d3;
            d_exact = d_exact3;
            Mesh = Mesh3;
        case 4
            d = d4;
            d_exact = d_exact4;
            Mesh = Mesh4;
        case 5
            d = d5;
            d_exact = d_exact5;
            Mesh = Mesh5;
        case 6
            d = d6;
            d_exact = d_exact6;
            Mesh = Mesh6;
    end
    
    % Mesh size
    switch Mesh.nsd
        case 1
            h(sim) = max(Mesh.coords)/Mesh.ne;
            Mesh.xdofs = Mesh.DOF;
            Mesh.ydofs = zeros(length(Mesh.DOF),1);
            Mesh.xdofs = Mesh.DOF;
            Mesh.ydofs = zeros(length(Mesh.DOF),1);
        case 2
            gcoords = Mesh.coords(Mesh.conn(1,:),:);
            h(sim) = sqrt(polyarea(gcoords(:,1),gcoords(:,2)));
    end
    
    % Calculate error norms
    eL2_num = 0;
    eL2_den = 0;
    
    % loop over elements
    for i = 1:Mesh.ne
        % element connectivity
        conn_e = Mesh.conn(i,:);
        % global coordinates
        gcoords = Mesh.coords(conn_e,:);
                
        % loop over integration points
        for ip = 1:Quad.nq
            % N matrices
            N = getN(Mesh, Quad, ip);
            % N derivatives
            dN = getdN(Mesh, Quad, ip);
            % Jacobian matrix
            J = dN*gcoords;
            % Jacobian determinant
            Jdet = det(J);
            
            % approximated displacement at quadrature point
            switch Mesh.nsd
                case 1
                    % approximated displacement at quadrature point
                    uxh = N*d(Mesh.xdofs(conn_e)');
                    uyh = zeros(length(uxh));
                    % exact displacement at quadrature point
                    uxe = N*d_exact(Mesh.xdofs(conn_e));
                    uye = zeros(length(uxe));

                case 2
                    % approximated displacement at quadrature point
                    uxh = N*d(Mesh.xdofs(conn_e)');
                    uyh = N*d(Mesh.ydofs(conn_e)');
                    % exact displacement at quadrature point
                    uxe = N*d_exact(Mesh.xdofs(conn_e)');
                    uye = N*d_exact(Mesh.ydofs(conn_e)');
            end

            % L2 norm
            eL2_num = eL2_num + [uxh - uxe , uyh - uye] * [uxh - uxe ; uyh - uye] * Quad.w(ip,1) * Jdet;
            eL2_den = eL2_den + [uxe , uye] * [uxe ; uye] * Quad.w(ip,1) * Jdet;
        end
    end
    
    eL2(sim) = sqrt(eL2_num);
end

% Determine slope of L2 norm
% pL2 = polyfit(log(h(1:2)), log(eL2p(1:2)),1);
% m_L2 = pL2(1);

figure;
loglog(h,eL2,'g-o', 'LineWidth', 1.5);
xlabel('Mesh size (m)');
ylabel('L2-norm');
title('Convergence for time');


end