function porousMedia_v1()
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
% E = 3e4; % elasticity modulus [Pa]
% k = 1e-10; % intrinsic permeability [m2]
% mu = 1e-3; % dynamic viscosity [Pa s]

E = 1;
k = 1;
mu = 1;

%% Mesh parameters
ne = 2; % number of elements
ndof_u_e = 3; % number of DOFs for u per element
ndof_p_e = 2; % number of DOFs for p per element
ndof_u = 5; % number of DOFs for u
ndof_p = 3; % number of DOFs for p
ndof_e = 5;  % number of DOFs per element
ndof = 8;   % total number of dofs

L = 1; % column size [m]
h = L/ne;   % element size [m]

coords = (0:h:L);

nq = 2;     % number of integration points

%% Solution parameters
t = 0;  % initial simulation time
dt = 1;  % time step
tend = 10;   % final simulation time
n = 1; % step counter
it = 0; % iteration counter
tol = 1e-3; % tolerance for the iterations

%% Initialize variables
u_old = zeros(ndof_u, 1);

fu = zeros(ndof_u, 1); % prescribed traction
fp = zeros(ndof_p, 1); % prescribed flux

plot_time = (0:dt:tend);
plot_u = zeros(1, length(plot_time));
plot_udot = zeros(1, length(plot_time));
plot_p = zeros(1, length(plot_time));

%% Boundary conditions
traction = -1e3; % applied traction [N]
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

%% Solve system

while t < tend
    % auxiliar terms
    fpbar = fp + (Kup.') * u_old /dt;
    
    % assemble full matrices
    A = [Kuu  -Kup;
        (Kup.')/dt  Kpp];
    b = [fu; fpbar];
    
    % solve system
    %    [Ldec,Udec,P] = lu(A);
    %    y = Ldec\(P*b);
    %    x = Udec\y;
    
       x = bicgstab(A,b);
    
    % ACCORDING TO MASTER'S NOTES
%     R = zeros(ndof,1); % residuals vector
%     x = zeros(ndof,1); % total displacements vector
%     
%     lambda = 0; % load factor
%     lambda_inc = 0.1; % load factor increment
%     
%     du = A\b; % displacements increment
%     Rnorm = 1000;
%     
%     while Rnorm > tol
%         if it == 0
%             dr = zeros(ndof,1); % residual increment for first iteration
%         else
%             dr = A\R; % residual increment for iterations > 1
%         end
%         
%         dx = lambda_inc*du + dr; % total variables increment
%         
%         % update variables
%         lambda = lambda + lambda_inc;
%         x = x + dx;
%         
%         % compute residuals
%         R = lambda*b - A*x;
%         
%         Rnorm = norm(R);
%         
%         it = it + 1; % update iteration counter
%     end
    
    
    % ACCORDING TO BRUCE'S CODE
    
%     x = A\b;
    Rnorm = 1000;
%     A11 = Kuu;
%     A12 = A(ndof_u+1:ndof, 1:ndof_u);
%     A21 = A(ndof_u+1:ndof, 1:ndof_u);
%     A22 = A(ndof_u+1:ndof, ndof_u+1:ndof);
%     
    
    while Rnorm > tol
        R = A*x - b; % residual
        
        Ru = R(1:ndof_u,1); % residual for displacement u
        Rp = R(ndof_u+1:ndof,1); % residual for pressure p
        
        u_norm = abs(traction); % value to normalize u
        p_norm = 1; % value to normalize p (would be zero, so use 1)
        
        Ru_norm = norm(R(1:ndof_u,1))*(1/u_norm); % norm of the residuals vector (u)
        Rp_norm = norm(R(ndof_u+1:ndof,1))*(1/p_norm); % norm of the residuals vector (p)
                
        Rnorm = sqrt(Ru_norm^2 + Rp_norm^2); % total residual norm
        
%         u_norm2 = mean(abs(diag(Kuu)));
%         p_norm2 = mean(abs(diag(Kpp)));
%         A11 = Kuu*(1/u_norm2);
%         A12 = (-Kup)*(1/u_norm2);
%         A21 = ((Kup.')/dt)*(1/p_norm2);
%         A22 = Kpp*(1/p_norm2);
%         A = [A11 A12; A21 A22];
%         
%         R(1:ndof_u,1) = R(1:ndof_u,1)*(1/u_norm2);
%         R(ndof_u+1:ndof,1) = R(ndof_u+1:ndof,1)*(1/p_norm2);
        
        dx = -A\R; % ????????????????????????? HERE WE HAVE A PROBLEM (Bruce's code: normalize all matrices for inversion

        dx(ndof_u, 1) = 0;    % no displacements at the bottom
        dx(ndof_u+1, 1) = 0; % drained top
        
        % update variables
        x = x + dx;
        
        % update iteration counter
        it = it + 1;
    end
    
    % store variables
    u = x(1 : ndof_u, 1);
    p = x(ndof_u+1 : ndof_u + ndof_p, 1);
    
    % enforce essential BCs
    u(ndof_u, 1) = 0;    % no displacements at the bottom
    p(1, 1) = 0; % drained top
    
    % enforce natural BCs
    fp(ndof_p, 1) = 0;    % impervious bottom
    
    % compute pressure variation in time
    udot = (u - u_old)./dt;
    
    plot_p(1,n) = p(ndof_p,1);
    
    % update variables
    u_old = u;
    t = t + dt;
    n = n + 1;
    
end

figure;
plot(plot_time, plot_p);

end