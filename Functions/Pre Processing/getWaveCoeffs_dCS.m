function [Ap_dCS, Bp_dCS, Cp_dCS, As_dCS, Bs_dCS, Cs_dCS] = getWaveCoeffs_dCS(Material, w)
% Compute polynomial coefficients for the wave solution for Biot (BT)
% theory
% ------------------------------------------------------------------------
% Based on equation from Zhao (2020): Effects of petrophysical parameters
% on attenuation and dispersion of seismic waves in the simplified
% poroelastic theory
% Return to de la Cruz and Spanos: include oscillations in porosity
% ------------------------------------------------------------------------

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

%% Auxiliar constants
aps = rhos-rho12/(1-eta0);
bps = rho12/(1-eta0);
cps = Ks*deltas/(1-eta0)-Ks-4*mus/3;
dps = Ks*deltaf/(1-eta0);
hps = muf*eta0^2/k/(1-eta0);

apf = rhof-rho12/eta0;
bpf = rho12/eta0;
cpf = Kf*deltaf/eta0-Kf;
dpf = Kf*deltas/eta0;
hpf = xif*deltaf/eta0-xif-4*muf/3;
mpf = xif*deltas/eta0;
rpf = muf*eta0/k;

ass = rhos-rho12/(1-eta0);
bss = rho12/(1-eta0);
css = muf*eta0^2/k/(1-eta0);
dss = mus;

asf = rhof-rho12/eta0;
bsf = rho12/eta0;
csf = muf*eta0/k;
dsf = muf;

%% Compressional P wave
Ap_dCS = cps*cpf - dps*dpf + (-cps*hpf + dps*mpf)*1i*w; 
Bp_dCS = (aps*cpf + hps*hpf + bps*dpf + dps*bpf - hps*mpf)*w.^2 +...
    (-aps*hpf - bps*mpf)*1i*w.^3 + ...
    (cps*apf + cps*rpf + hps*cpf - dps*rpf - hps*dpf)*1i*w;
Cp_dCS = -bps*bpf*w.^4 - hps*apf*w.^2 + ...
    + (aps*apf + aps*rpf + bps*rpf + hps*bpf)*1i*w.^3;

%% Shear S wave
As_dCS = -dss*dsf*1i*w;
Bs_dCS = -(css*dsf+dss*asf)*w.^2 + ass*dsf*1i*w.^3 - dss*csf*1i*w;
Cs_dCS = (ass*asf - bss*bsf)*w.^4 + (ass*csf + css*asf + bss*csf + css*bsf)*1i*w.^3;

end