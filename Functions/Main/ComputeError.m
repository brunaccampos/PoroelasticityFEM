function ComputeError()
% ------------------------------------------------------------------------
% Compute error using the Le norm of the displacements
% ------------------------------------------------------------------------

clearvars
clear, clc
close all

%% Number of sims
nsims = 4;

%% Mesh 1
% load file
load Results_m1_space_Boone.mat
% mesh
MeshP1 = MeshP;
MeshU1 = MeshU;
% FEM approximation
dp1 = p;
du1 = u;
% exact solution
dp_exact1 = p_an;
du_exact1 = u_an;

%% Mesh 2
load Results_m2_space_Boone.mat
% mesh
MeshP2 = MeshP;
MeshU2 = MeshU;
% FEM approximation
dp2 = p;
du2 = u;
% exact solution
dp_exact2 = p_an;
du_exact2 = u_an;

%% Mesh 3
load Results_m3_space_Boone.mat
% mesh
MeshP3 = MeshP;
MeshU3 = MeshU;
% FEM approximation
dp3 = p;
du3 = u;
% exact solution
dp_exact3 = p_an;
du_exact3 = u_an;

%% Mesh 4
load Results_m4_space_Boone.mat
% mesh
MeshP4 = MeshP;
MeshU4 = MeshU;
% FEM approximation
dp4 = p;
du4 = u;
% exact solution
dp_exact4 = p_an;
du_exact4 = u_an;

%% Mesh 5
% load Results_m5_space_Boone.mat
% % mesh
% MeshP5 = MeshP;
% MeshU5 = MeshU;
% % FEM approximation
% dp5 = p;
% du5 = u;
% % exact solution
% dp_exact5 = p_an;
% du_exact5 = u_an;




% shear modulus [GPa]
Material.G = 6;
% Poisson's ratio
Material.nu = 0.2;
% elasticity modulus [GPa]
Material.E = 2 * Material.G * (1 + Material.nu);
% porous media permeability [m2/GPa s]
Material.kf = 2e-2;
% dynamic viscosity [GPa s]
Material.mu = 1e-12;
% intrinsic permeability [m2]
Material.k = Material.kf * Material.mu;
% Biot's coefficient
Material.alpha = 1;
% fluid bulk modulus [GPa]
Material.Kf = 3;
% solid bulk modulus [GPa]
Material.Ks = 36;
% fluid bulk viscosity [GPa s]
Material.xif = 2.8e-12; % (Quiroga-Goode, 2005)
% material porosity
Material.n = 0.19;
% 1/Q (related to storage coefficient)
Material.Minv = (Material.alpha - Material.n)/Material.Ks + Material.n/Material.Kf;

% additional coefficients for analytical result
% Lame constant [GPa]
Material.lambda = Material.E * Material.nu/((1+Material.nu)*(1-2*Material.nu));
% gravitational acceleration [m/s2]
Material.g = 9.81;
% fluid density [10^9 kg/m3]
Material.rho_f = 1000e-9;
% hydraulic conductivity [m/s]
Material.kh = Material.kf * Material.rho_f * Material.g;





% column length
L = max(MeshU.coords);

% coordinates vector
% xu = MeshU.coords;
% xp = MeshP.coords;

syms x

% traction
P = 1e-6;
% time
t = 10;

% material parameters
mu = Material.G;
M = 1/Material.Minv;
alpha = Material.alpha;
Ku = Material.lambda + 2*mu/3 + alpha^2*M;
cm = 1/(Material.lambda + 2*mu);
c = Material.kh/(Material.rho_f*Material.g*(Material.Minv + alpha^2*cm));

% initial displacement after instantaneous load
u0 = P*(L-x)/(Ku + 4*mu/3);
% initial pressure after instantaneous load
p0 = alpha*M*P/(Ku + 4*mu/3);

% number of terms for Fourier expansion
N = 100;

% auxiliar terms
auxp = 0;
auxu = 0;

% loop over N
for m = 0:N
   auxp = auxp + (1/(2*m+1)) * exp(-(2*m+1)^2*pi()^2*c*t/(4*L^2)) * sin((2*m+1)*pi().*x/(2*L));
   auxu = auxu + (1/(2*m+1)^2) * exp(-(2*m+1)^2*pi()^2*c*t/(4*L^2)) * cos((2*m+1)*pi().*x/(2*L));
end

p_an = 4*p0.*auxp/pi();
u_an = cm*p0*(L-x -8*L.*auxu/pi()^2) + u0;





%% Compute L2 error norm
% initialize variables
eL2p = zeros(nsims,1);
eL2u = zeros(nsims,1);
h = zeros(nsims,1);
normu = zeros(nsims,1);
normp = zeros(nsims,1);

% loop over meshes
for sim = 1:nsims
    
    switch sim
        case 1
            dp = dp1;
            du = du1;
            dp_exact = dp_exact1;
            du_exact = du_exact1;
            MeshP = MeshP1;
            MeshU = MeshU1;
            normu(sim,1) = norm(du-du_exact,2);
            normp(sim,1) = norm(dp-dp_exact,2);
        case 2
            dp = dp2;
            du = du2;
            dp_exact = dp_exact2;
            du_exact = du_exact2;
            MeshP = MeshP2;
            MeshU = MeshU2;
            normu(sim,1) = norm(du-du_exact,2);
            normp(sim,1) = norm(dp-dp_exact,2);
        case 3
            dp = dp3;
            du = du3;
            dp_exact = dp_exact3;
            du_exact = du_exact3;
            MeshP = MeshP3;
            MeshU = MeshU3;
            normu(sim,1) = norm(du-du_exact,2);
            normp(sim,1) = norm(dp-dp_exact,2);
        case 4
            dp = dp4;
            du = du4;
            dp_exact = dp_exact4;
            du_exact = du_exact4;
            MeshP = MeshP4;
            MeshU = MeshU4;
            normu(sim,1) = norm(du-du_exact,2);
            normp(sim,1) = norm(dp-dp_exact,2);
        case 5
            dp = dp5;
            du = du5;
            dp_exact = dp_exact5;
            du_exact = du_exact5;
            MeshP = MeshP5;
            MeshU = MeshU5;
            normu(sim,1) = norm(du-du_exact,2);
            normp(sim,1) = norm(dp-dp_exact,2);
        case 6
            dp = dp6;
            du = du6;
            dp_exact = dp_exact6;
            du_exact = dp_exact6;
            MeshP = MeshP6;
            MeshU = MeshU6;
            normu(sim,1) = norm(du-du_exact,2);
            normp(sim,1) = norm(dp-dp_exact,2);
    end
    
    % Mesh size
    switch MeshP.nsd
        case 1
            h(sim) = max(MeshP.coords)/MeshP.ne;
            MeshP.xdofs = MeshP.DOF;
            MeshP.ydofs = zeros(length(MeshP.DOF),1);
            MeshU.xdofs = MeshU.DOF;
            MeshU.ydofs = zeros(length(MeshU.DOF),1);
        case 2
            gcoordsp = MeshP.coords(MeshP.conn(1,:),:);
            h(sim) = sqrt(polyarea(gcoordsp(:,1),gcoordsp(:,2)));
    end
    
    % exact solution
%     d_exact = @(x) x.^5 - x.^4;

    Control.nqU = 8;
    MeshU.field = 'u';
    Quad = GlobalQuad(MeshU, Control);
    
    % Calculate error norms
    eL2p_num = 0;
    eL2p_den = 0;
    eL2u_num = 0;
    eL2u_den = 0;
    
    % loop over elements
    for i = 1:MeshP.ne
        % element connectivity
        connp_e = MeshP.conn(i,:);
        connu_e = MeshU.conn(i,:);
        % global coordinates
        gcoordsp = MeshP.coords(connp_e,:);
                
        % loop over integration points
        for ip = 1:Quad.nq
            % N matrices
            Np = getN(MeshP, Quad, ip);
            Nu = getN(MeshU, Quad, ip);
            % N derivatives
            dNp = getdN(MeshP, Quad, ip);
            % Jacobian matrix
            J = dNp*gcoordsp;
            % Jacobian determinant
            Jdet = det(J);
            
            % approximated displacement at quadrature point
            switch MeshP.nsd
                case 1
                    % approximated displacement at quadrature point
                    uxh_p = Np*dp(MeshP.xdofs(connp_e)');
                    uyh_p = zeros(length(uxh_p));
                    uxh_u = Nu*du(MeshU.xdofs(connu_e)');
                    uyh_u = zeros(length(uxh_u));
                    % exact displacement at quadrature point
                    uxe_p = Np*dp_exact(MeshP.xdofs(connp_e));
                    uye_p = zeros(length(uxe_p));
                    uxe_u = Nu*du_exact(MeshU.xdofs(connu_e));
                    uye_u = zeros(length(uxe_u));

%                     uxe_p = eval(subs(p_an, Np*MeshP.coords(connp_e)));
%                     uxe_u = eval(subs(u_an, Nu*MeshU.coords(connu_e)));

%                     uxe_u = eval(subs(d_exact, Nu*MeshU.coords(connu_e)));

                case 2
                    % approximated displacement at quadrature point
                    uxh_p = Np*dp(MeshP.xdofs(connp_e)');
                    uyh_p = Np*dp(MeshP.ydofs(connp_e)');
                    uxh_u = Nu*du(MeshU.xdofs(connu_e)');
                    uyh_u = Nu*du(MeshU.ydofs(connu_e)');
                    % exact displacement at quadrature point
                    uxe_p = Np*dp_exact(MeshP.xdofs(connp_e)');
                    uye_p = Np*dp_exact(MeshP.ydofs(connp_e)');
                    uxe_u = Nu*du_exact(MeshU.xdofs(connu_e)');
                    uye_u = Nu*du_exact(MeshU.ydofs(connu_e)');
            end

            % L2 norm pressure
            eL2p_num = eL2p_num + [uxh_p - uxe_p , uyh_p - uye_p] * [uxh_p - uxe_p ; uyh_p - uye_p] * Quad.w(ip,1) * Jdet;
            eL2p_den = eL2p_den + [uxe_p , uye_p] * [uxe_p ; uye_p] * Quad.w(ip,1) * Jdet;
            % L2 norm displacement
            eL2u_num = eL2u_num + [uxh_u - uxe_u , uyh_u - uye_u] * [uxh_u - uxe_u ; uyh_u - uye_u] * Quad.w(ip,1) * Jdet;
            eL2u_den = eL2u_den + [uxe_u , uye_u] * [uxe_u ; uye_u] * Quad.w(ip,1) * Jdet;
        end
    end
    
    eL2p(sim) = sqrt(eL2p_num/eL2p_den);
    eL2u(sim) = sqrt(eL2u_num/eL2u_den);
end

% Determine slope of L2 norm
pL2p = polyfit(log(h(1:2)), log(eL2p(1:2)),1);
m_L2p = pL2p(1);

pL2u = polyfit(log(h(1:2)), log(eL2u(1:2)),1);
m_L2u = pL2u(1);

pL2p2 = polyfit(log(h), log(normp),1);
m_L2p2 = pL2p2(1);

pL2u2 = polyfit(log(h), log(normu),1);
m_L2u2 = pL2u2(1);


%% Step 2 - Calculate the slope of each curve

figure;
loglog(h,eL2p,'m-o', 'LineWidth', 1.5);
% hold on
% loglog(h,normp,'b-o', 'LineWidth', 1.5);
% hold off
xlabel('Mesh size (m)');
ylabel('L2-norm');
title(sprintf('Convergence for pressure order %.2f', m_L2p));

figure;
loglog(h,eL2u,'g-o', 'LineWidth', 1.5);
% hold on
% loglog(h,normu,'k-o', 'LineWidth', 1.5);
% hold off
xlabel('Mesh size (m)');
ylabel('L2-norm');
title(sprintf('Convergence for displacement order %.2f', m_L2u));


end