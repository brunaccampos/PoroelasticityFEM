% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function [p_an, u_an] = getAnSol_coupledIncomp_Symb(Control, Material, MeshU, ~, BC)
% Analytical solution for column consolidation
% Valid for incompressible solid and fluid materials

% Reference: Korsawe, J., Starke, G., Wang, W., & Kolditz, O. (2006). 
% Finite element analysis of poro-elastic consolidation in porous media: 
% Standard and mixed approaches. Computer Methods in Applied Mechanics and 
% Engineering, 195(9-12), 1096-1115.

% column length
L = max(MeshU.coords);

% coordinates vector
syms x

% additional material parameters
mu = Material.M(1).E/(2*(1+Material.M(1).nu)); % shear modulus [Pa]
lambda = Material.M(1).E*Material.M(1).nu/((1+Material.M(1).nu)*(1-2*Material.M(1).nu)); % Lamé constant [Pa]

aux1 = 0; % u
aux2 = 0; % sigma
aux3 = 0; % p

% number of terms considered for the Fourier expansion
N = 100;

% dimensionless time
td = (lambda + 2*mu)*Material.M(1).k * Control.t / (Material.M(1).muf * L^2);

% loop of N terms
for n = 0:N
    M = pi()*(2*n+1)/2;
    aux1 = aux1 + 2* cos(M*x)* exp(-M^2 * td)/M^2;
    aux2 = aux2 + 2* sin(M*x)* exp(-M^2 * td)/M;
    aux3 = aux3 + 2* sin(M*x)* exp(-M^2 * td)/M;
end

% displacement
ud = 1 - x - aux1;
u_an = ud * BC.pointLoadValue * L/(lambda + 2*mu);

% pressure
pd = aux3;
p_an = pd * BC.pointLoadValue;

end