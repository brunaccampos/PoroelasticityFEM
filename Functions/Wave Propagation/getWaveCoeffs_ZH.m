% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function [Ap_ZH, Bp_ZH, Cp_ZH, As_ZH, Bs_ZH, Cs_ZH] = getWaveCoeffs_ZH(Material, w)
% Compute polynomial coefficients for the wave solution for Biot (BT)
% theory
% Based on equation from Zhao (2020): Effects of petrophysical parameters
% on attenuation and dispersion of seismic waves in the simplified
% poroelastic theory

rhof = Material.rhof; % fluid density
muf = Material.muf; % fluid dynamic viscosity
Kf = Material.Kf; % fluid bulk modulus
xif = Material.xif; % fluid bulk viscosity
rhos = Material.rhos; % solid density 
mus = Material.mus; % solid shear modulus
Ks = Material.Ks; % solid bulk modulus
rho12 = Material.rho12; % coupled density
eta0 = Material.eta0; % porosity
k = Material.k; % intrinsic permeability

%% Auxiliar constants
vps2 = (Ks+4*mus/3)/rhos;
aps = muf*eta0^2/(k*rhos*(1-eta0));
cps = rho12/(1-eta0)/rhos;

c2 = Kf/rhof;
bpf = muf*eta0/k/rhof;
cpf = rho12/eta0/rhof;
mpf = (xif+4*muf/3)/rhof;

css = mus/rhos;
dss = muf*eta0^2/(k*rhos*(1-eta0));

dsf = muf*eta0/k/rhof;
csf = muf/rhof;

%% Compressional P wave
Ap_ZH = c2*vps2-mpf*vps2*1i*w;
Bp_ZH = -(c2*aps+bpf*vps2)*1i*w+(-vps2-c2+c2*cps-mpf*aps+cpf*vps2)*w.^2+(mpf-mpf*cps)*1i*w.^3;
Cp_ZH = (1-cps-cpf)*w.^4+(aps+bpf)*1i*w.^3;

%% Shear S wave
As_ZH = -csf*css*1i*w;
Bs_ZH = (-css-csf*dss+css*cpf)*w.^2+(csf-csf*cps)*1i*w.^3-dsf*css*1i*w;
Cs_ZH = (1-cps-cpf)*w.^4+(dss+dsf)*1i*w.^3;

end