clear
clc
clear vars
close all

%% Material parameters
n = 0.19; % porosity [-]
k = 1.88e-13; % intrinsic permeability [m2]
Kf = 3300e6; % fluid bulk modulus [Pa]
Ks = 36000e6; % solid bulk modulus [Pa]
mu = 1e-3; % dynamic viscosity [Pa s]
rhos = 2600; % solid density [kg/m^3]
rhof = 1000; % fluid density [kg/m^3]

%% Spanos parameters
deltaF = n*(1-n)/(Ks*(n/Kf + (1-n)/Ks));
deltaS = n*(1-n)/(Kf*(n/Kf + (1-n)/Ks));

%% Biot parameters
A = (1-n-deltaS) * Ks - (1-n)*mu*2/3; % Lamé constant
N = (1-n) * mu; % Lamé constant
Q = Kf * deltaS; % coupling term
R = n * Kf * (1-deltaF/n);
P = A + 2*N;
H = P + R + 2*Q;
rho11 = (1-n)*rhos;
rho22 = n*rhof;
rho = rho11 + rho22;

% nondimensional parameters
sigma11 = P/H;
sigma22 = R/H;
sigma12 = Q/H;
gamma11 = rho11/rho;
gamma22 = rho22/rho;
gamma12 = 0;

% characteristic frequency
fc = mu*n^2/(k*2*pi*rho*(gamma12+gamma22));
f_min = 20e3;
f_max = 2e6;
f = (f_min:f_max);

a = sigma11*sigma22-sigma12^2;
b = -(sigma22*gamma11 + sigma11*gamma22 - 2*sigma12*gamma12);
c = (gamma11*gamma22 - gamma12^2);

p = [a b c];
z = sort(roots(p));
zeta1 = z(1,1) - 1;
zeta2 = z(2,1) - 1;

% frequency range
% f = (0:0.001:0.15);
f_range = f/fc;

% phase velocity - rotational wave
vR = 1+ (1/8) *(4*gamma22 - (gamma12+gamma22)^2) * f_range.^2;
% attenuation coefficient - rotational wave
attR = (1/2)*(gamma12 + gamma22) * f_range.^2;

% phase velocity - dilatational waves 1st kind
vD1 = 1 - 0.5*(f_range.^2) * ((sigma11*sigma22 - sigma12)^2 / (gamma12+gamma22)^2) * zeta1 * zeta2 * (zeta1+zeta2+0.5*zeta1*zeta2);
% attenuation - dilatational waves 1st kind
attD1 = 0.5*abs(zeta1*zeta2) * (sigma11*sigma22 - sigma12^2)/(gamma12 + gamma22) * (f_range.^2);

% phase velocity - dilatational waves 2nd kind
vD2 = sqrt(2*f_range * (sigma11*sigma22 - sigma12^2)/(gamma12 + gamma22));
% attenuation - dilatational waves 2nd kind
attD2 = sqrt(0.5*f_range * (gamma12 + gamma22)/(sigma11 * sigma22 - sigma12^2));

%% Plots
% velocity - rotational
figure;
subplot(3,2,1);
plot(f_range, vR, 'b', 'LineWidth', 2);
hold on
xlabel('Frequency range');
ylabel('Velocity');
title('Rotational waves - Phase velocity');
hold off

% attenuation - rotational
subplot(3,2,2);
plot(f_range, attR, 'm', 'LineWidth', 2);
hold on
xlabel('Frequency range');
ylabel('Attenuation');
title('Rotational waves - Attenuation');
hold off

% velocity - dilatational 1st
subplot(3,2,3);
plot(f_range, vD1, 'b', 'LineWidth', 2);
hold on
xlabel('Frequency range');
ylabel('Velocity');
title('Dilatational waves I - Phase velocity');
hold off

% attenuation - dilatational 1st
subplot(3,2,4);
plot(f_range, attD1, 'm', 'LineWidth', 2);
hold on
xlabel('Frequency range');
ylabel('Attenuation');
title('Dilatational waves I - Attenuation');
hold off

% velocity - dilatational 2nd
subplot(3,2,5);
plot(f_range, vD2, 'b', 'LineWidth', 2);
hold on
xlabel('Frequency range');
ylabel('Velocity');
title('Dilatational waves II - Phase velocity');
hold off

% attenuation - dilatational 2nd
subplot(3,2,6);
plot(f_range, attD2, 'm', 'LineWidth', 2);
hold on
xlabel('Frequency range');
ylabel('Attenuation');
title('Dilatational waves II - Attenuation');
hold off