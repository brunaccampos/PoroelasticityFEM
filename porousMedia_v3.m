function porousMedia_v3()
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
mu = 1e-3; % dynamic viscosity [Pa s]

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

xd_u = coords_u./L;
xd_p = coords_p./L;

top_node_u = find(coords_u == min(coords_u));
bottom_node_u = find(coords_u == max(coords_u));

top_node_p = find(coords_p == min(coords_p));
bottom_node_p = find(coords_p == max(coords_p));

fixed_u = [bottom_node_u];
fixed_p = [top_node_p];
fixed_fp = [bottom_node_p];

dof_u = (1:ndof_u);
dof_p = (1:ndof_p);

free_u = setdiff(dof_u, fixed_u);
free_p = setdiff(dof_p, fixed_p);
free_fp = setdiff(dof_p, fixed_fp);

%% Solution parameters
t = 0;  % initial simulation time
dt = 1;  % time step
tend = 10;   % final simulation time
n = 1; % step counter

%% Initialize variables
u_old = zeros(ndof_u, 1);

fu = zeros(ndof_u, 1); % prescribed traction
fp = zeros(ndof_p, 1); % prescribed flux

u = zeros(ndof_u,1);
p = zeros(ndof_p,1);

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
        Kpp_e = Kpp_e + (k/mu) * (Bp.') * Bp * J * w_csi(ip,1);

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


Kpu = (Kup.')/dt;

Kuu_EE = Kuu(fixed_u, fixed_u);
Kuu_EF = Kuu(fixed_u,free_u);
Kuu_FE = Kuu(free_u, fixed_u);
Kuu_FF = Kuu(free_u, free_u);

Kup_EE = Kup(fixed_u, fixed_p);
Kup_EF = Kup(fixed_u,free_p);
Kup_FE = Kup(free_u, fixed_p);
Kup_FF = Kup(free_u, free_p);

Kpu_EE = Kpu(fixed_p, fixed_u);
Kpu_EF = Kpu(fixed_p,free_u);
Kpu_FE = Kpu(free_p, fixed_u);
Kpu_FF = Kpu(free_p, free_u);

Kpp_EE = Kpp(fixed_p, fixed_p);
Kpp_EF = Kpp(fixed_p,free_p);
Kpp_FE = Kpp(free_p, fixed_p);
Kpp_FF = Kpp(free_p, free_p);

% reassemble
KEE = [Kuu_EE, -Kup_EE;
    Kpu_EE, Kpp_EE];
KEF = [Kuu_EF, -Kup_EF;
    Kpu_EF, Kpp_EF];
KFE = [Kuu_FE, -Kup_FE;
    Kpu_FE, Kpp_FE];
KFF = [Kuu_FF, -Kup_FF;
    Kpu_FF, Kpp_FF];

fuE = fu(free_u);
fuF = fu(fixed_u);

%% Solve system

while t < tend
    % analytical solution
    [~, p_an, u_an] = getAnalyticResult (traction, E, nu, mu, k, xd_u, xd_p, L, t);

    % auxiliar terms
    fpbar = fp + (Kup.') * u_old /dt;
    %     fubar = fu*t/tend; % incremental loading

    fpE = fpbar(free_p);
    fpF = fpbar(fixed_p);

    uE = u(fixed_u);
    pE = p(fixed_p);

    dE = [uE;pE];
    fE = [fuE; fpE];

    % solve for unknown variables
    dF = KFF\(fE - KFE *dE);

    % solve for unknown reactions
    fF = KEE*dE + KEF*dF;

    % store variables
    u = [dF(1:length(free_u),1); uE];
    p = [pE; dF(length(free_u)+1 : length(free_u) + length(free_p),1)];

    % update variables
    u_old = u;
    t = t + dt;
    n = n + 1;

end

% displacement plot
figure;
plot(coords_u,u,'b','LineWidth',2);
hold on
plot(coords_u, u_an,'m--','LineWidth',2);
xlabel('Column depth (m)');
ylabel('u (m)');
title(sprintf('Displacement at t = %.0f s', tend));
legend('Numerical', 'Analytical');
hold off
exportgraphics(gcf,'Displacement.png','Resolution',300)

% pressure plot
figure;
plot(coords_p,p,'b','LineWidth',2);
hold on
plot(coords_p,p_an,'m--','LineWidth',2);
xlabel('Column depth (m)');
ylabel('p (Pa)');
title(sprintf('Pressure at t = %.0f s', tend));
legend('Numerical', 'Analytical');
hold off
exportgraphics(gcf,'Pressure.png','Resolution',300)

end