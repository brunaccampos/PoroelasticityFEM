% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function [Ap_dCS, Bp_dCS, Cp_dCS, As_dCS, Bs_dCS, Cs_dCS] = getWaveCoeffs_dCS(Material, w)
% Compute polynomial coefficients for the wave solution for Biot (BT)
% theory
% Based on equation from Zhao (2020): Effects of petrophysical parameters
% on attenuation and dispersion of seismic waves in the simplified
% poroelastic theory
% Return to de la Cruz and Spanos: include oscillations in porosity

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

%% Auxiliar constants
vps2 = (Ks+4*mus*(1+c)/3)/rhos;
aps = muf*eta0^2/(k*rhos*(1-eta0));
cps = rho12/(1-eta0)/rhos;
dps = Ks*deltas/rhos/(1-eta0);
hps = Ks*deltaf/rhos/(1-eta0);

c2 = Kf/rhof;
bpf = muf*eta0/k/rhof;
cpf = rho12/eta0/rhof;
dpf = Kf*deltas/rhof/eta0;
hpf = Kf*deltaf/rhof/eta0;
mpf = xif*deltaf/rhof/eta0-xif/rhof-4*muf/3/rhof; % new in dCS
rpf = xif*deltas/rhof/eta0 + 4*(1-eta0)*muf*c/3/rhof/eta0; % new in dCS

css = mus/rhos;
dss = muf*eta0^2/(k*rhos*(1-eta0));

dsf = muf*eta0/k/rhof;
csf = muf/rhof; % new in dCS

%% Compressional P wave
Ap_dCS = c2*vps2-vps2*hpf-dps*c2+dps*hpf-hps*dpf + (vps2*mpf-dps*mpf+hps*rpf)*1i*w;
Bp_dCS = (-c2+hpf-vps2+vps2*cpf+dps-dps*cpf+cps*c2-cps*hpf+hps*cpf+cps*dpf+aps*mpf-aps*rpf)*w.^2 + ...
    (-vps2*bpf+dps*bpf-aps*c2+aps*hpf-hps*bpf-aps*dpf)*1i*w + ...
    (cps*mpf-cps*rpf-mpf)*1i*w.^3;
Cp_dCS = (1-cps-cpf)*w.^4+1i*(aps+bpf)*w.^3;

%% Shear S wave
As_dCS = -css*csf*1i*w;
Bs_dCS = (-css+css*cpf-dss*csf)*w.^2-dsf*css*1i*w+(csf-cps*csf)*1i*w.^3;
Cs_dCS = (1-cps-cpf)*w.^4+(dss+dsf)*1i*w.^3;

end