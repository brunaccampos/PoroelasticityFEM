function [p_an, u_an] = getAnalyticResult_v3(MeshU, MeshP, Control)
% ------------------------------------------------------------------------
% Analytical solution for column consolidation (Terzaghi problem)
% Reference: Ferronato (2010)
% ------------------------------------------------------------------------
% version 3: same as v2 but with symbolic function. Data from 1D Boone
% ------------------------------------------------------------------------

%% Material properties - Boone (1990)
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

% symbolic coordinate
syms x

% traction
P = 1e-6;
% time
t = Control.t;

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

end