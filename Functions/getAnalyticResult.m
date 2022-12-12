% Analytical solution for column consolidation

function [sigma_an, p_an, u_an] = getAnalyticResult (Material, Mesh, BC, Control)

% dimensionless coordinates
xd_u = Mesh.coords_u./Mesh.L;
xd_p = Mesh.coords_p./Mesh.L;

% additional material parameters
G = Material.E/(2*(1+Material.nu)); % shear modulus [Pa]
lambda = Material.E*Material.nu/((1+Material.nu)*(1-2*Material.nu)); % Lamé constant [Pa]

aux1 = 0; % u
aux2 = 0; % sigma
aux3 = 0; % p

% number of terms considered for the Fourier expansion
N = 100;

% dimensionless time
td = (lambda + 2*G)*Material.k * Control.t / (Material.mu * Mesh.L^2);

% loop of N terms
for n = 0:N
    M = pi()*(2*n+1)/2;
    aux1 = aux1 + 2* cos(M*xd_u)* exp(-M^2 * td)/M^2;
    aux2 = aux2 + 2* sin(M*xd_u)* exp(-M^2 * td)/M;
    aux3 = aux3 + 2* sin(M*xd_p)* exp(-M^2 * td)/M;
end

% displacement
ud = 1 - xd_u - aux1;
u_an = ud*BC.traction*Mesh.L/(lambda + 2*G);

% stress
sigmad = -1 + aux2;
sigma_an = sigmad*BC.traction;

% pressure
pd = aux3;
p_an = pd*BC.traction;

end