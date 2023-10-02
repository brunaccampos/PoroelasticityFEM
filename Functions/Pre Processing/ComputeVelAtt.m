function ComputeVelAtt()
% Compute velocity and attenuation for compressional and shear waves using
% different porous media models
% ------------------------------------------------------------------------
% Models to compare:
% Biot (1941, 1956a): porosity considered implicitly, no fluid viscous
% dissipation terms
% Spanos (2002): explicit porosity equation, includes fluid viscous
% dissipation terms
% Zhao (2020): simplified version from Spanos where oscillations in
% porosity are disregarded
% ------------------------------------------------------------------------

%% Material parameters - fixed
Material.rhof = 1050; % fluid density [kg/m3]
Material.muf = 1e-3; % fluid dynamic viscosity [Pa s]
Material.Kf = 2.2e9; % fluid bulk modulus [Pa]
Material.xif = 2.8e-3; % fluid bulk viscosity [Pa s]

Material.rhos = 2650; % solid density [kg/m3]
Material.mus = 23e9; % solid shear modulus [Pa]
Material.Ks = 33e9; % solid bulk modulus [Pa]

Material.rho12 = -83; % coupled density [kg/m3]

%% Material parameters - study variation
Material.eta0 = 0.15; % porosity [-]
Material.k = 10*9.8692331e-16; % permeability [m2] Note: 1D = 1e-12 m2, 1mD = 1e-15 m2

%% Frequency array
w = 1e2:100:1e7; % frequency [Hz]

%% Compute polynomial constants
% Biot (BT) theory
[Ap_BT, Bp_BT, Cp_BT, As_BT, Bs_BT, Cs_BT] = getWaveCoeffs_BT(Material);
% de la Cruz and Spanos (dCS) theory
[Ap_dCS, Bp_dCS, Cp_dCS, As_dCS, Bs_dCS, Cs_dCS] = getWaveCoeffs_dCS(Material);
% Zhao - simplified de la Cruz and Spanos - (ZH) theory
[Ap_ZH, Bp_ZH, Cp_ZH, As_ZH, Bs_ZH, Cs_ZH] = getWaveCoeffs_ZH(Material);
% zero vector
Z = zeros(length(w),1);

% regroup constants - P wave
Ap = [Ap_BT; Ap_dCS; Ap_ZH];
Bp = [Bp_BT; Bp_dCS; Bp_ZH];
Cp = [Cp_BT; Cp_dCS; Cp_ZH];
% regroup constants - S wave
As = [As_BT; As_dCS; As_ZH];
Bs = [Bs_BT; Bs_dCS; Bs_ZH];
Cs = [Cs_BT; Cs_dCS; Cs_ZH];

%% Compute polynomial roots
% loop over models
for m = 1:3
    % polynomial 4th order P wave
    pol_p = [Ap(m)' Z Bp' Z Cp'];
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
    
    
end


%% Compressional P wave
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