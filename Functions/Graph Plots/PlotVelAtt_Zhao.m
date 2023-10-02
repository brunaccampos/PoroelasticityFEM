function PlotVelAtt_Zhao()
% Plot phase velocity and attenuation for compressional and shear waves
% ------------------------------------------------------------------------
% Based on equation from Zhao (2020): Effects of petrophysical parameters
% on attenuation and dispersion of seismic waves in the simplified
% poroelastic theory
% ------------------------------------------------------------------------

clear vars
clc
close all

%% Material parameters - fixed
rhof = 1050; % fluid density [kg/m3]
muf = 1e-3; % fluid dynamic viscosity [Pa s]
Kf = 2.2e9; % fluid bulk modulus [Pa]
xif = 2.8e-3; % fluid bulk viscosity [Pa s]

rhos = 2650; % solid density [kg/m3]
mus = 23e9; % solid shear modulus [Pa]
Ks = 33e9; % solid bulk modulus [Pa]

rho12 = -83; % coupled density [kg/m3]

%% Material parameters - study variation
eta0 = 0.15; % porosity [-]
k = 10*9.8692331e-16; % permeability [m2] Note: 1D = 1e-12 m2, 1mD = 1e-15 m2

%% Frequency array
w = 1e2:100:1e7; % frequency [Hz]

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
% coefficients
Ap = c2*vps2-1i*apf*vps2*w;
Bp = -1i*(c2*aps+bpf*vps2)*w-(vps2+c2-c2*cps+apf*aps-cpf*vps2)*w.^2+1i*(apf-apf*cps)*w.^3;
Cp = (1-cps-cpf)*w.^4+1i*(aps+bpf)*w.^3;
Z = zeros(length(w),1);
% polynomial 4th order
pol_p = [Ap' Z Bp' Z Cp'];
% initialize matrices
roots_p = zeros(length(w),4);
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
As = -1i*csf*css*w;
Bs = -(css+csf*dss-css*cpf)*w.^2+1i*(csf-csf*cps)*w.^3-1i*dsf*css*w;
Cs = (1-cps-cpf)*w.^4+1i*(dss+dsf)*w.^3;
Z = zeros(length(w),1);
% polynomial 4th order
pol_s = [As' Z Bs' Z Cs'];
% initialize matrices
roots_s = zeros(length(w),4);
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
semilogx(w,vp(:,2),'b--', 'LineWidth', 1.5);
hold on
grid on
xlabel('Frequency [Hz]');
ylabel('Velocity [m/s]');
title('First P wave');
hold off

% second P wave
nexttile;
semilogx(w,vp(:,1), 'g--', 'LineWidth', 1.5);
hold on
grid on
xlabel('Frequency [Hz]');
ylabel('Velocity [m/s]');
title('Second P wave');
hold off

% first S wave
nexttile;
semilogx(w,vs(:,2), 'k--', 'LineWidth', 1.5);
hold on
grid on
xlabel('Frequency [Hz]');
ylabel('Velocity [m/s]');
title('First S wave');
hold off

% second S wave
nexttile;
semilogx(w,vs(:,1), 'r--', 'LineWidth', 1.5);
hold on
grid on
xlabel('Frequency [Hz]');
ylabel('Velocity [m/s]');
title('Second S wave');
hold off

%% Plots attenuation
% first P wave
nexttile;
semilogx(w,attp(:,2),'b', 'LineWidth', 1.5);
hold on
grid on
xlabel('Frequency [Hz]');
ylabel('Attenuation');
title('First P wave');
hold off

% second P wave
nexttile;
semilogx(w,attp(:,1), 'g', 'LineWidth', 1.5);
hold on
grid on
xlabel('Frequency [Hz]');
ylabel('Attenuation');
title('Second P wave');
hold off

% first S wave
nexttile;
semilogx(w,atts(:,2), 'k', 'LineWidth', 1.5);
hold on
grid on
xlabel('Frequency [Hz]');
ylabel('Attenuation');
title('First S wave');
hold off

% second S wave
nexttile;
semilogx(w,atts(:,1), 'r', 'LineWidth', 1.5);
hold on
grid on
xlabel('Frequency [Hz]');
ylabel('Attenuation');
title('Second S wave');
hold off

end