function [p_an, u_an] = getAnSol_coupledComp(Control, Material, MeshU, MeshP, BC)
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
xu = MeshU.coords;
xp = MeshP.coords;

% traction
P = BC.pointLoadValue;
% time
t = Control.t;

% material parameters
mu = Material.mu;
M = 1/Material.Minv;
alpha = Material.alpha;
Ku = Material.lambda + 2*mu/3 + alpha^2*M;
cm = 1/(Material.lambda + 2*mu);
c = Material.kh/(Material.rhof*Material.g*(Material.Minv + alpha^2*cm));

% initial displacement after instantaneous load
u0 = P*(L-xu)/(Ku + 4*mu/3);
% initial pressure after instantaneous load
p0 = alpha*M*P/(Ku + 4*mu/3);

% number of terms for Fourier expansion
N = 1000;

% auxiliar terms
auxp = 0;
auxu = 0;

% loop over N
for m = 0:N
   auxp = auxp + (1/(2*m+1)) * exp(-(2*m+1)^2*pi()^2*c*t/(4*L^2)) * sin((2*m+1)*pi().*xp/(2*L));
   auxu = auxu + (1/(2*m+1)^2) * exp(-(2*m+1)^2*pi()^2*c*t/(4*L^2)) * cos((2*m+1)*pi().*xu/(2*L));
end

p_an = 4*p0.*auxp/pi();
u_an = cm*p0*(L-xu -8*L.*auxu/pi()^2) + u0;

end