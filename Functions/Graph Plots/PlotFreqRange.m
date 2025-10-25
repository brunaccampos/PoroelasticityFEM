% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function PlotFreqRange(Material, f0, range)
% Plot frequency range for given material parameters

%% Range to plot
switch range
    case 1 % seismic
        f_min = 1;
        f_max = 100;
    case 2 % acoustic
        f_min = 20;
        f_max = 20e3;
    case 3 % ultrasonic
        f_min = 20e3;
        f_max = 2e6;
end
% frequency vector
f = (f_min:10:f_max);

%% Material parameters
n = Material.eta0; % porosity [-]
k = Material.k; % intrinsic permeability [m2]
Kf = Material.Kf; % fluid bulk modulus [Pa]
Ks = Material.Ks; % solid bulk modulus [Pa]
muf = Material.muf; % dynamic viscosity [Pa s]
rhos = Material.rhos; % solid density [kg/m^3]
rhof = Material.rhof; % fluid density [kg/m^3]

%% Spanos parameters
deltaf = Material.deltaf;
deltas = Material.deltas;

%% Biot parameters
A = (1-n-deltas) * Ks - (1-n)*muf*2/3; % Lamé constant
N = (1-n) * muf; % Lamé constant
Q = Kf * deltas; % coupling term
R = n * Kf * (1-deltaf/n);
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
fc = muf*n^2/(k*2*pi*rho*(gamma12+gamma22));

% current frequency ratio
ratio = f0/fc;
fprintf('Current frequency: %.2f Hz \n', f0);
fprintf('Characteristic frequency: %.2f Hz \n', fc);
fprintf('Ratio f/fc = %.2f \n', ratio);

a = sigma11*sigma22-sigma12^2;
b = -(sigma22*gamma11 + sigma11*gamma22 - 2*sigma12*gamma12);
c = (gamma11*gamma22 - gamma12^2);

p = [a b c];
z = sort(roots(p));
zeta1 = z(1,1) - 1;
zeta2 = z(2,1) - 1;

% frequency range
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
subplot(2,3,1);
loglog(f_range, vR, 'b', 'LineWidth', 2);
grid on
hold on
xline(ratio, '--', 'LineWidth', 1.5);
xlabel('Frequency range f/f_c [-]');
ylabel('Velocity [-]');
title('Rotational waves - Phase velocity');
hold off

% attenuation - rotational
subplot(2,3,4);
loglog(f_range, attR, 'm', 'LineWidth', 2);
grid on
hold on
xline(ratio, '--', 'LineWidth', 1.5);
xlabel('Frequency range f/f_c [-]');
ylabel('Attenuation [-]');
title('Rotational waves - Attenuation');
hold off

% velocity - dilatational 1st
subplot(2,3,2);
loglog(f_range, vD1, 'b', 'LineWidth', 2);
grid on
hold on
xline(ratio, '--', 'LineWidth', 1.5);
xlabel('Frequency range f/f_c [-]');
ylabel('Velocity [-]');
title('Dilatational waves I - Phase velocity');
hold off

% attenuation - dilatational 1st
subplot(2,3,5);
loglog(f_range, attD1, 'm', 'LineWidth', 2);
grid on
hold on
xline(ratio, '--', 'LineWidth', 1.5);
xlabel('Frequency range f/f_c [-]');
ylabel('Attenuation [-]');
title('Dilatational waves I - Attenuation');
hold off

% velocity - dilatational 2nd
subplot(2,3,3);
loglog(f_range, vD2, 'b', 'LineWidth', 2);
grid on
hold on
xline(ratio, '--', 'LineWidth', 1.5);
xlabel('Frequency range f/f_c [-]');
ylabel('Velocity [-]');
title('Dilatational waves II - Phase velocity');
hold off

% attenuation - dilatational 2nd
subplot(2,3,6);
loglog(f_range, attD2, 'm', 'LineWidth', 2);
grid on
hold on
xline(ratio, '--', 'LineWidth', 1.5);
xlabel('Frequency range f/f_c [-]');
ylabel('Attenuation [-]');
title('Dilatational waves II - Attenuation');
hold off

end