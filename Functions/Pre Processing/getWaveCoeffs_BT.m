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
deltas = Material.deltas; % compliance coefficient solid
deltaf = Material.deltaf; % compliance coefficient fluid
eta0 = Material.eta0; % porosity
k = Material.k; % intrinsic permeability

%% Auxiliar constants
vps2 = (Ks+4*mus/3)/rhos;
aps = muf*eta0^2/(k*rhos*(1-eta0));
cps = rho12/(1-eta0)/rhos;
dps = Ks*deltas/rhos/(1-eta0);
hps = Ks*deltaf/rhos/(1-eta0);

c2 = Kf/rhof;
bpf = muf*eta0/k/rhof;
cpf = rho12/eta0/rhof;
dpf = Kf*deltas/rhof/eta0;
hpf = Kf*deltaf/rhof/eta0;

css = mus/rhos;
dss = muf*eta0^2/(k*rhos*(1-eta0));
dsf = muf*eta0/k/rhof;


%% Compressional P wave
Ap_BT = (c2*vps2-vps2*hpf-dps*c2+dps*hpf-hps*dpf)*ones(1,length(w));
Bp_BT = (-c2*aps-bpf*vps2+dps*bpf+aps*hpf-hps*bpf-aps*dpf)*1i*w +...
    (-vps2-c2+c2*cps+cpf*vps2+hpf+dps-dps*cpf-cps*hpf+hps*cpf+cps*dpf)*w.^2;
Cp_BT = (1-cps-cpf)*w.^4+1i*(aps+bpf)*w.^3;

%% Shear S wave
As_BT = zeros(1,length(w));
Bs_BT = -(css-css*cpf)*w.^2-1i*dsf*css*w;
Cs_BT = (1-cps-cpf)*w.^4+1i*(dss+dsf)*w.^3;

end