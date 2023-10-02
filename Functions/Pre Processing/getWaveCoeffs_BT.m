function [Ap_BT, Bp_BT, Cp_BT, As_BT, Bs_BT, Cs_BT] = getWaveCoeffs_BT(Material, w)
% Compute polynomial coefficients for the wave solution for Biot (BT)
% theory
% ------------------------------------------------------------------------
% Based on equation from Zhao (2020): Effects of petrophysical parameters
% on attenuation and dispersion of seismic waves in the simplified
% poroelastic theory
% Return to Biot: assume xif = 0 (fluid bulk viscosity)
% ------------------------------------------------------------------------

rhof = Material.rhof; % fluid density
muf = Material.muf; % fluid dynamic viscosity
Kf = Material.Kf; % fluid bulk modulus
xif = 0; % fluid bulk viscosity
rhos = Material.rhos; % solid density 
mus = Material.mus; % solid shear modulus
Ks = Material.Ks; % solid bulk modulus
rho12 = Material.rho12; % coupled density
eta0 = Material.eta0; % porosity
k = Material.k; % intrinsic permeability

%% Auxiliar constants
c2 = Kf/rhof;
vps2 = (Ks+4*mus/3)/rhos;
aps = muf*eta0^2/(k*rhos*(1-eta0));
apf = (xif+4*muf/3)/rhof;
bpf = muf*eta0/k/rhof;
cps = rho12/(1-eta0)/rhos;
cpf = rho12/eta0/rhof;
css = mus/rhos;
csf = muf/rhof;
dss = muf*eta0^2/(k*rhos*(1-eta0));
dsf = muf*eta0/k/rhof;

%% Compressional P wave
Ap_BT = c2*vps2-1i*apf*vps2*w;
Bp_BT = -1i*(c2*aps+bpf*vps2)*w-(vps2+c2-c2*cps+apf*aps-cpf*vps2)*w.^2+1i*(apf-apf*cps)*w.^3;
Cp_BT = (1-cps-cpf)*w.^4+1i*(aps+bpf)*w.^3;

%% Shear S wave
As_BT = -1i*csf*css*w;
Bs_BT = -(css+csf*dss-css*cpf)*w.^2+1i*(csf-csf*cps)*w.^3-1i*dsf*css*w;
Cs_BT = (1-cps-cpf)*w.^4+1i*(dss+dsf)*w.^3;

end