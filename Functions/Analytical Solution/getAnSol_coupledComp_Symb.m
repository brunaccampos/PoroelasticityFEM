function [p_an, u_an] = getAnSol_coupledComp_Symb(Control, Material, MeshU, ~, BC)
% ------------------------------------------------------------------------
% Analytical solution for column consolidation (Terzaghi problem)
% Valid for compressible solid and fluid materials
% ------------------------------------------------------------------------
% Reference: Ferronato, M., Castelletto, N., & Gambolati, G. (2010). A 
% fully coupled 3-D mixed finite element model of Biot consolidation. 
% Journal of Computational Physics, 229(12), 4813-4830.
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
mu = Material.M(1).mu;
M = 1/Material.M(1).Minv;
alpha = Material.M(1).alpha;
Ku = Material.M(1).lambda + 2*mu/3 + alpha^2*M;
cm = 1/(Material.M(1).lambda + 2*mu);
c = Material.M(1).kh/(Material.M(1).rhof*Material.M(1).g*(Material.M(1).Minv + alpha^2*cm));

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