function ComputeVelAtt()
% Compute velocity and attenuation for compressional and shear waves using
% different porous media models
% ------------------------------------------------------------------------
% Models to compare:
% 1- Biot (1941, 1956a): porosity considered implicitly, no fluid viscous
%   dissipation terms
% 2- Zhao (2020): simplified version from Spanos where oscillations in
%   porosity are disregarded
% 3- Spanos (2002): explicit porosity equation, includes fluid viscous
%   dissipation terms
% ------------------------------------------------------------------------

clear vars
clc
close all

%% Material parameters - fixed
Material.rhof = 1050; % fluid density [kg/m3]
Material.muf = 1e-3; % fluid dynamic viscosity [Pa s]
Material.Kf = 2.2e9; % fluid bulk modulus [Pa]
Material.xif = 2.8e-3; % fluid bulk viscosity [Pa s]

Material.rhos = 2650; % solid density [kg/m3]
Material.mus = 23e9; % solid shear modulus [Pa]
Material.Ks = 33e9; % solid bulk modulus [Pa]

Material.rho12 = -83; % coupled density [kg/m3]
Material.alpha = 0.79; % Biot coefficient [-]

%% Material parameters - study variation
Material.eta0 = 0.15; % porosity [-]
Material.k = 10*9.8692331e-16; % permeability [m2] Note: 1D = 1e-12 m2, 1mD = 1e-15 m2

%% Material parameters - dCS model
% porosity effective pressure coefficient (Spanos, 1989)
% Material.n = 0; % lower limit
Material.n = 1; % return to Biot
% Material.n = Material.Ks/Material.Kf; % upper limit

Minv = Material.eta0/Material.Kf + (Material.alpha-Material.eta0)/Material.Ks;
Mstarinv = Minv-(1-Material.n)*(Material.alpha-Material.eta0)/Material.Ks;
Mstar = 1/Mstarinv;

% porosity equation coefficients
Material.deltas = (Material.alpha-Material.eta0)*Material.eta0*Mstar/Material.Kf;
Material.deltaf = (Material.alpha-Material.eta0)*Material.eta0*Mstar*Material.n/Material.Ks;

%% Frequency array
w = 1e2:1e3:1e8; % frequency [Hz]

%% Compute polynomial constants
% Biot (BT) theory
[Ap_BT, Bp_BT, Cp_BT, As_BT, Bs_BT, Cs_BT] = getWaveCoeffs_BT(Material, w);
% Zhao - simplified de la Cruz and Spanos - (ZH) theory
[Ap_ZH, Bp_ZH, Cp_ZH, As_ZH, Bs_ZH, Cs_ZH] = getWaveCoeffs_ZH(Material, w);
% de la Cruz and Spanos (dCS) theory
[Ap_dCS, Bp_dCS, Cp_dCS, As_dCS, Bs_dCS, Cs_dCS] = getWaveCoeffs_dCS(Material, w);
% zero vector
Z = zeros(length(w),1);

% regroup constants - P wave
Ap = [Ap_BT; Ap_ZH; Ap_dCS];
Bp = [Bp_BT; Bp_ZH; Bp_dCS];
Cp = [Cp_BT; Cp_ZH; Cp_dCS];
% regroup constants - S wave
As = [As_BT; As_ZH; As_dCS];
Bs = [Bs_BT; Bs_ZH; Bs_dCS];
Cs = [Cs_BT; Cs_ZH; Cs_dCS];

%% Initialize variables
vp = zeros(length(w),6);
vs = zeros(length(w),6);
attp = zeros(length(w),6);
atts = zeros(length(w),6);

%% Compute polynomial roots
% loop over models
for m = 1:3
    % polynomial 4th order P wave
    pol_p = [Ap(m,:)' Z Bp(m,:)' Z Cp(m,:)'];
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
    vp(:, 2*m-1:2*m) = w'./real(kp);
    % attenuations
    attp(:, 2*m-1:2*m) = abs(2*imag(kp)./real(kp));
    
    % polynomial 4th order S wave
    pol_s = [As(m,:)' Z Bs(m,:)' Z Cs(m,:)'];
    % initialize matrices
    roots_s = zeros(length(w),4);
    ks = zeros(length(w),2);
    % loop over frequency values
    for i = 1:length(w)
        % roots
        aux = roots(pol_s(i,:));
        roots_s(i,1:length(aux)) = aux;
        % physical solutions
        ks(i,1:length(aux)/2) = roots_s(i, real(roots_s(i,:))>0);
    end
    % phase velocities
    vs(:, 2*m-1:2*m) = w'./real(ks);
    % attenuations
    atts(:, 2*m-1:2*m) = abs(2*imag(ks)./real(ks));   
end

% initialize figure
figure;
tiledlayout(2,4);

%% Plots phase velocity
% first P wave
nexttile;
semilogx(w,vp(:,2),'k-', 'LineWidth', 1.5); % BT
hold on
grid on
xlabel('Frequency [Hz]');
ylabel('Velocity [m/s]');
title('First P wave');
semilogx(w,vp(:,4),'b--', 'LineWidth', 1.5); % ZH
semilogx(w,vp(:,6),'r:', 'LineWidth', 1.5); % dCS
legend('BT', 'ZH', 'dCS');
hold off

% second P wave
nexttile;
semilogx(w,vp(:,1), 'k-', 'LineWidth', 1.5); % BT
hold on
grid on
xlabel('Frequency [Hz]');
ylabel('Velocity [m/s]');
title('Second P wave');
semilogx(w,vp(:,3), 'b--', 'LineWidth', 1.5); % ZH
semilogx(w,vp(:,5), 'r:', 'LineWidth', 1.5); % dCS
legend('BT', 'ZH', 'dCS');
hold off

% first S wave
nexttile;
semilogx(w,vs(:,1), 'k-', 'LineWidth', 1.5); % BT
hold on
grid on
xlabel('Frequency [Hz]');
ylabel('Velocity [m/s]');
title('First S wave');
semilogx(w,vs(:,4), 'b--', 'LineWidth', 1.5); % ZH
semilogx(w,vs(:,6), 'r:', 'LineWidth', 1.5); % dCS
legend('BT', 'ZH', 'dCS');
hold off

% second S wave
nexttile;
semilogx(w,vs(:,3), 'b--', 'LineWidth', 1.5); % ZH
hold on
grid on
xlabel('Frequency [Hz]');
ylabel('Velocity [m/s]');
title('Second S wave');
semilogx(w,vs(:,5), 'r:', 'LineWidth', 1.5); % dCS
legend('ZH', 'dCS');
hold off

%% Plots attenuation
% first P wave
nexttile;
semilogx(w,attp(:,2),'k-', 'LineWidth', 1.5); % BT
hold on
grid on
xlabel('Frequency [Hz]');
ylabel('Attenuation');
title('First P wave');
semilogx(w,attp(:,4),'b--', 'LineWidth', 1.5); % ZH
semilogx(w,attp(:,6),'r:', 'LineWidth', 1.5); % dCS
legend('BT', 'ZH', 'dCS');
hold off

% second P wave
nexttile;
semilogx(w,attp(:,1), 'k-', 'LineWidth', 1.5); % BT
hold on
grid on
xlabel('Frequency [Hz]');
ylabel('Attenuation');
title('Second P wave');
semilogx(w,attp(:,3), 'b--', 'LineWidth', 1.5); % ZH
semilogx(w,attp(:,5), 'r:', 'LineWidth', 1.5); % dCS
legend('BT', 'ZH', 'dCS');
hold off

% first S wave
nexttile;
semilogx(w,atts(:,1), 'k-', 'LineWidth', 1.5); % BT
hold on
grid on
xlabel('Frequency [Hz]');
ylabel('Attenuation');
title('First S wave');
semilogx(w,atts(:,4), 'b--', 'LineWidth', 1.5); % ZH
semilogx(w,atts(:,6), 'r:', 'LineWidth', 1.5); % dCS
legend('BT', 'ZH', 'dCS');
hold off

% second S wave
nexttile;
semilogx(w,atts(:,3), 'b--', 'LineWidth', 1.5); % ZH
hold on
grid on
xlabel('Frequency [Hz]');
ylabel('Attenuation');
title('Second S wave');
semilogx(w,atts(:,5), 'r:', 'LineWidth', 1.5); % dCS
legend('ZH', 'dCS');
hold off

end