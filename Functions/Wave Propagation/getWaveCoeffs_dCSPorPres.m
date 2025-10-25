% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function [AdCS, BdCS, CdCS] = getWaveCoeffs_dCSPorPres(Material, w)
% Compute polynomial coefficients for the wave solution for Biot (BT)
% theory

%% Material parameters
rhof = Material.rhof; % fluid density
muf = Material.muf; % fluid dynamic viscosity
Kf = Material.Kf; % fluid bulk modulus
xif = Material.xif; % fluid bulk viscosity
rhos = Material.rhos; % solid density 
mus = Material.mus; % solid shear modulus
Ks = Material.Ks; % solid bulk modulus
rho12 = Material.rho12; % coupled density
deltas = Material.deltas; % compliance coefficient solid
deltaf = Material.deltaf; % compliance coefficient fluid
eta0 = Material.eta0; % porosity
k = Material.k; % intrinsic permeability
c = Material.c; % micro heterogeneity
mum = (1-eta0)*(1+c)*mus; % macroscopic shear modulus

%% Auxiliar constants
A = -rhof/eta0 + rho12/eta0^2 + rho12*(1/deltas-deltaf/deltas/eta0)/eta0;
B = muf*eta0/k*(-1/eta0-1/deltas+deltaf/deltas/eta0);
C = 4*muf/3/eta0 - (1-eta0)*muf/eta0*(mum/(1-eta0)/mus-1)*4/3*(1/deltas-deltaf/deltas/eta0);
D = rhof/Kf - rho12/eta0/Kf + rho12*deltaf/eta0/deltas/Kf;
E = - muf*eta0/k*(-1/Kf+deltaf/deltas*Kf);
F = -(xif+4*muf/3)/Kf - (1-eta0)*muf/eta0*(mum/(1-eta0)/mus-1)*4/3*deltaf/deltas/Kf;
G = rhos*(1/deltas-deltaf/deltas/eta0) - rho12/(1-eta0)/eta0 - rho12/(1-eta0)*(1/deltas-deltaf/deltas/eta0);
H = muf*eta0^2/k/(1-eta0)*(1/eta0+1/deltas-deltaf/deltas/eta0);
I = -Ks*(1/deltas-deltaf/deltas/eta0) + Ks/(1-eta0) - 4*mum/3/(1-eta0)*(1/deltas-deltaf/deltas/eta0);
J = rhos*deltaf/deltas/Kf + rho12/(1-eta0)/Kf - rho12*deltaf/(1-eta0)/deltas/Kf;
K = muf*eta0^2/k/(1-eta0)*(-1/Kf+deltaf/deltas/Kf);
L = -Ks*deltaf/deltas/Kf - 4*mum/3/(1-eta0)*deltaf/deltas/Kf;

%% Compressional P waves
AdCS = C*L*1i*w - F*1i*w*I + I;
BdCS = C*J*1i*w.^3 + C*K*1i^2*w.^2 - F*G*1i*w.^3 - F*H*1i^2*w.^2 - A*L*w.^2 - ...
    B*L*1i*w + D*w.^2*I + E*1i*w*I + G*w.^2 + H*1i*w;
CdCS = -A*J*w.^4 - A*K*1i*w.^3 - B*J*1i*w.^3 - B*K*1i^2*w.^2 + D*G*w.^4 + ...
    D*H*1i*w.^3 + E*G*1i*w.^3 + E*H*1i^2*w.^2;

end