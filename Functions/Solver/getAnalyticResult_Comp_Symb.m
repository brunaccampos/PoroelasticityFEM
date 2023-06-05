function [p_an, u_an] = getAnalyticResult_Comp_Symb(Material, MeshU, BC, Control)
% ------------------------------------------------------------------------
% Analytical solution for column consolidation (Terzaghi problem)
% Valid for compressible solid and fluid materials
% ------------------------------------------------------------------------
% Reference: Ferronato (2010)
% ------------------------------------------------------------------------

% column length
L = max(MeshU.coords);

% coordinates vector
syms x

% traction
P = BC.pointLoadValue;
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