function porousMedia_v2()
% Porous media simulation
% Bruna Campos
% September 2021

% Based on Korsawe (2006) model

% Assumptions/conventions:
% - stress is positive for tension
% - boundary condition for force is based on total stress
% - no acceleration terms for solid or fluid
% - solid velocity is neglected
% - fluid and solid grains are incompressible, but porous medium itself is
% compressible

clearvars
clear
clc

%% Material parameters
E = 3e4; % elasticity modulus [Pa]
nu = 0.2; % Poisson's ratio
k = 1e-10; % intrinsic permeability [m2]
eta = 1e-3; % dynamic viscosity [Pa s]

% sanity check
% E = 1;
% nu = 1;
% k = 1;
% eta = 1;

%% Mesh parameters (change)
ne = 10; % number of elements
L = 1; % column size [m]

%% Mesh (do not change)
ndof_u_e = 3; % number of DOFs for u per element
ndof_p_e = 2; % number of DOFs for p per element

ndof_u = ne * (ndof_u_e - 1) + 1; % number of DOFs for u
ndof_p = ne * (ndof_p_e - 1) + 1; % number of DOFs for p

h = L/ne;   % element size [m]

coords_u = (0:(L/(ndof_u-1)):L); % nodal coordinates for displacement field
coords_p = (0:(L/(ndof_p-1)):L); % nodal coordinates for pressure field

nq = 2;     % number of integration points

%% Solution parameters
t = 0;  % initial simulation time
dt = 1;  % time step
tend = 10;   % final simulation time
n = 1; % step counter

%% Initialize variables
u_old = zeros(ndof_u, 1);

fu = zeros(ndof_u, 1); % prescribed traction
fp = zeros(ndof_p, 1); % prescribed flux

%% Boundary conditions
traction = 1e3; % applied traction [N]
fu(1,1) = traction; % Neumann BC

%% Integration points
csi = zeros(nq,1);
w_csi = zeros(nq,1);

% coordinates
csi(1,1) = -0.5773502692;
csi(2,1) = 0.5773502692;

% weights
w_csi(:,1) = 1;

% Jacobian
J = h/2;

%% Assemble matrices

% initialize global matrices
Kuu = zeros(ndof_u, ndof_u);
Kup = zeros(ndof_u,ndof_p);
Kpp = zeros (ndof_p, ndof_p);

% loop over elements
for e = 1:ne
    % initialize local matrices
    Kuu_e = zeros(ndof_u_e, ndof_u_e);
    Kup_e = zeros(ndof_u_e,ndof_p_e);
    Kpp_e = zeros (ndof_p_e, ndof_p_e);

    % loop over integration points
    for ip = 1:nq
        % shape functions for u
        Nu = zeros(1, 3);
        Nu(1,1) = 0.5*(csi(ip,1)^2 - csi(ip,1));
        Nu(1,2) = 1 - csi(ip,1)^2;
        Nu(1,3) = 0.5*(csi(ip,1)^2 + csi(ip,1));

        % derivative shape functions for u
        dNu = zeros(1,3);
        dNu(1,1) = csi(ip,1) - 1/2;
        dNu(1,2) = -2*csi(ip,1);
        dNu(1,3) = csi(ip,1) + 1/2;

        % B matrix for u
        Bu = dNu./J;

        % shape functions for p
        Np = zeros(1, 2);
        Np(1,1) = 0.5*(1 - csi(ip,1));
        Np(1,2) = 0.5*(1 + csi(ip,1));

        % derivative shape functions for p
        dNp = zeros(1,2);
        dNp(1,1) = -1/2;
        dNp(1,2) = 1/2;

        % B matrix for p
        Bp = dNp./J;

        % assemble local matrices
        Kuu_e = Kuu_e + (Bu.') * E * Bu * J * w_csi(ip,1);
        Kup_e = Kup_e + (Bu.') * Np * J * w_csi(ip,1);
        Kpp_e = Kpp_e + (k/eta) * (Bp.') * Bp * J * w_csi(ip,1);

    end

    % connectivity matrix for u
    Connu = zeros(ndof_u_e, ndof_u);
    for i = 1:ndof_u_e
        for j = 1:ndof_u
            if j == (ndof_u_e-1)*(e-1)+i
                Connu(i,j) = 1;
            end
        end
    end

    % connectivity matrix for p
    Connp = zeros(ndof_p_e, ndof_p);
    for i = 1:ndof_p_e
        for j = 1:ndof_p
            if j == (ndof_p_e-1)*(e-1)+i
                Connp(i,j) = 1;
            end
        end
    end

    % assemble global matrices
    Kuu = Kuu + (Connu.') * Kuu_e * Connu;
    Kup = Kup + (Connu.') * Kup_e * Connp;
    Kpp = Kpp + (Connp.') * Kpp_e * Connp;

end

%% Solve system

while t < tend
    % analytical solution
    xd = coords_u./L;
    [~, p_an, u_an] = getAnalyticResult (traction, E, nu, eta, k, xd, L, t);

    % auxiliar terms
    fpbar = fp + (Kup.') * u_old /dt;
    %     fubar = fu*t/tend; % incremental loading

    % assemble full matrices
    A = [Kuu  -Kup;
        (Kup.')/dt  Kpp];
    b = [fu; fpbar];

    % convert to sparse matrices
    As = sparse(A);
    bs = sparse(b);

    % incomplete LU factorization
    options.type = 'ilutp';
    options.droptol = 1; % ignore nondiagonal values smaller than this tolerance

    [Ldec,Udec] = ilu(As,options);

    % solve system
    x = bicgstab(As,bs,[],[], Ldec, Udec);

    % store variables
    u = x(1 : ndof_u, 1);
    p = x(ndof_u+1 : ndof_u + ndof_p, 1);

    % enforce essential BCs
    u(ndof_u, 1) = 0;    % no displacements at the bottom
    p(1, 1) = 0; % drained top

    % enforce natural BCs
    fp(ndof_p, 1) = 0;    % impervious bottom

    % update variables
    u_old = u;
    t = t + dt;
    n = n + 1;

end

% displacement plot
figure;
plot(coords_u,u);
hold on
plot(coords_u, u_an);
xlabel('Column depth (m)');
ylabel('u (m)');
title(sprintf('Displacement at t = %.0f s', tend));
legend('Numerical', 'Analytical');
hold off

% pressure plot
figure;
plot(coords_p,p);
hold on
plot(coords_u,p_an);
xlabel('Column depth (m)');
ylabel('p (Pa)');
title(sprintf('Pressure at t = %.0f s', tend));
legend('Numerical', 'Analytical');
hold off

end