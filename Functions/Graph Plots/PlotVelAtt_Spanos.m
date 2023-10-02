function PlotVelAtt_Spanos()
% Plot phase velocity and attenuation for compressional and shear waves
% ------------------------------------------------------------------------
% Based on equation from Zhao (2020): Effects of petrophysical parameters
% on attenuation and dispersion of seismic waves in the simplified
% poroelastic theory
% Adapted to account for porosity variation
% ------------------------------------------------------------------------

clear vars
clc
close all

%% Material parameters - fixed
rhof = 1050; % fluid density [kg/m3]
muf = 1e-3; % fluid dynamic viscosity [Pa s]
xif = 2.8e-3; % fluid bulk viscosity [Pa s]
Kf = 2.2e9; % fluid bulk modulus [Pa]
rhos = 2650; % solid density [kg/m3]
mus = 23e9; % solid shear modulus [Pa]
Ks = 33e9; % solid bulk modulus [Pa]
rho12 = -83; % coupled density [kg/m3]

%% Material parameters - study variation
eta0 = 0.15; % porosity [-]
k = 10*1e-15; % permeability [m2] Note: 1D = 1e-12 m2, 1mD = 1e-15 m2

%% Material parameters - dCS model
alpha = 0.79; % Biot coefficient [-]

% porosity effective pressure coefficient (Spanos, 1989)
% n = 0; % lower limit
n = 1; % return to Biot
% n = Material.Ks/Material.Kf; % upper limit

Minv = eta0/Kf+(alpha-eta0)/Ks;
Mstarinv = Minv-(1-n)*(alpha-eta0)/Ks;
Mstar = 1/Mstarinv;

% porosity equation coefficients
deltas = (alpha-eta0)*eta0*Mstar/Kf;
deltaf = (alpha-eta0)*eta0*Mstar*n/Ks;

%% Frequency array
w = 1e2:100:1e7; % frequency [Hz]

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
% coefficients
Ap = cps*cpf - dps*dpf + (-cps*hpf + dps*mpf)*1i*w; 
Bp = (aps*cpf + hps*hpf + bps*dpf + dps*bpf - hps*mpf)*w.^2 +...
    (-aps*hpf - bps*mpf)*1i*w.^3 + ...
    (cps*apf + cps*rpf + hps*cpf - dps*rpf - hps*dpf)*1i*w;
Cp = -bps*bpf*w.^4 - hps*apf*w.^2 + ...
    + (aps*apf + aps*rpf + bps*rpf + hps*bpf)*1i*w.^3;
Z = zeros(length(w),1);
% polynomial 4th order
pol_p = [Ap' Z Bp' Z Cp'];
% initialize matrices
roots_p = zeros(length(w), 4);
kp = zeros(length(w),2);
% loop over frequency values
for i = 1:length(w)
    % roots
    roots_p(i,:) = roots(pol_p(i,:));
    % physical solutions
    kp(i,:) = roots_p(i, real(roots_p(i,:))>0);
end
% phase velocities
vp = w'./real(kp);
% attenuations
attp = abs(2*imag(kp)./real(kp));

%% Shear S wave
As = -dss*dsf*1i*w;
Bs = -(css*dsf+dss*asf)*w.^2 + ass*dsf*1i*w.^3 - dss*csf*1i*w;
Cs = (ass*asf - bss*bsf)*w.^4 + (ass*csf + css*asf + bss*csf + css*bsf)*1i*w.^3;
Z = zeros(length(w),1);
% polynomial 4th order
pol_s = [As' Z Bs' Z Cs'];
% initialize matrices
roots_s = zeros(length(w), 4);
ks = zeros(length(w),2);
% loop over frequency values
for i = 1:length(w)
    % roots
    roots_s(i,:) = roots(pol_s(i,:));
    % physical solutions
    ks(i,:) = roots_s(i, real(roots_s(i,:))>0);
end
% phase velocities
vs = w'./real(ks);
% attenuations
atts = abs(2*imag(ks)./real(ks));

% initialize figure
figure;
tiledlayout(2,4);

%% Plots phase velocity
% first P wave
nexttile;
semilogx(w,vp(:,1),'b--', 'LineWidth', 1.5);
hold on
grid on
xlabel('Frequency [Hz]');
ylabel('Velocity [m/s]');
title('First P wave');
hold off

% second P wave
nexttile;
semilogx(w,vp(:,2), 'g--', 'LineWidth', 1.5);
hold on
grid on
xlabel('Frequency [Hz]');
ylabel('Velocity [m/s]');
title('Second P wave');
hold off

% first S wave
nexttile;
semilogx(w,vs(:,1), 'k--', 'LineWidth', 1.5);
hold on
grid on
xlabel('Frequency [Hz]');
ylabel('Velocity [m/s]');
title('First S wave');
hold off

% second S wave
nexttile;
semilogx(w,vs(:,2), 'r--', 'LineWidth', 1.5);
hold on
grid on
xlabel('Frequency [Hz]');
ylabel('Velocity [m/s]');
title('Second S wave');
hold off

%% Plots attenuation
% first P wave
nexttile;
semilogx(w,attp(:,1),'b', 'LineWidth', 1.5);
hold on
grid on
xlabel('Frequency [Hz]');
ylabel('Attenuation');
title('First P wave');
hold off

% second P wave
nexttile;
semilogx(w,attp(:,2), 'g', 'LineWidth', 1.5);
hold on
grid on
xlabel('Frequency [Hz]');
ylabel('Attenuation');
title('Second P wave');
hold off

% first S wave
nexttile;
semilogx(w,atts(:,1), 'k', 'LineWidth', 1.5);
hold on
grid on
xlabel('Frequency [Hz]');
ylabel('Attenuation');
title('First S wave');
hold off

% second S wave
nexttile;
semilogx(w,atts(:,2), 'r', 'LineWidth', 1.5);
hold on
grid on
xlabel('Frequency [Hz]');
ylabel('Attenuation');
title('Second S wave');
hold off

end