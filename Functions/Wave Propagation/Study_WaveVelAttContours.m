% Study of wave velocity and attenuation: computation for different
% frequency values and range of material parameters
% ------------------------------------------------------------------------
% Models to compare:
% 1- Biot (1941, 1956a): porosity considered implicitly, no fluid viscous
%   dissipation terms
% 2- Zhao (2020): simplified version from Spanos where oscillations in
%   porosity are disregarded
% 3- Spanos (2002): explicit porosity equation, includes fluid viscous
%   dissipation terms
% ------------------------------------------------------------------------

% clear vars
% clc
% close all

%% Frequency array
w = logspace(0,6); % frequency [Hz]

%% Material parameters (Zhao 2020)
% Material.rhof = 1050; % fluid density [kg/m3]
% Material.muf = 1e-3; % fluid dynamic viscosity [Pa s]
% Material.Kf = 2.2e9; % fluid bulk modulus [Pa]
% Material.xif = 2.8e-3; % fluid bulk viscosity [Pa s]
%
% Material.rhos = 2650; % solid density [kg/m3]
% Material.mus = 23e9; % solid shear modulus [Pa]
% Material.Ks = 33e9; % solid bulk modulus [Pa]
%
% Material.rho12 = -83; % coupled density [kg/m3]
% Material.alpha = 0.79; % Biot coefficient [-]
% Material.eta0 = 0.15; % porosity [-]
% Material.k = 10*9.8692331e-16; % permeability [m2] Note: 1D = 1e-12 m2, 1mD = 1e-15 m2

%% Material parameters (Detournay 1993)
Material.rhof = 1000; % fluid density [kg/m3]
Material.muf = 1e-3; % fluid dynamic viscosity [Pa s]
Material.Kf = 3.3e9; % fluid bulk modulus [Pa]
Material.xif = 2.8e-3; % fluid bulk viscosity [Pa s]

Material.rhos = 2600; % solid density [kg/m3]
Material.mus = 27e9; % solid shear modulus [Pa]
Material.Ks = 36e9; % solid bulk modulus [Pa]

Material.rho12 = 0; % coupled density [kg/m3]
Material.alpha = 0.79; % Biot coefficient [-]
Material.eta0 = 0.19; % porosity [-]
Material.k = 1.88e-13; % permeability [m2] Note: 1D = 1e-12 m2, 1mD = 1e-15 m2

%% Material parameters (Spanos 1985)
% Material.rhof = 750; % fluid density [kg/m3]
% Material.muf = 1e4; % fluid dynamic viscosity [Pa s]
% Material.Kf = 3.34e10*9807e-9; % fluid bulk modulus [Pa] Note: nT/m2 = 9807e-9 Pa
% Material.xif = 0; % fluid bulk viscosity [Pa s]
%
% Material.rhos = 2650; % solid density [kg/m3]
% Material.mus = 1.5e10*9807e-9; % solid shear modulus [Pa] Note: nT/m2 = 9807e-9 Pa
% Material.Ks = 3; % solid bulk modulus [Pa]
%
% Material.rho12 = 0; % coupled density [kg/m3]
% Material.eta0 = 0.30; % porosity [-]
% Material.k = 1e-8; % permeability [m2] Note: 1D = 1e-12 m2, 1mD = 1e-15 m2
%
% % converting parameters
% R = 8.96e8*9807e-9;
% Q = 4.008e7*9807e-9;
% M = R/Material.eta0^2;
%
% Material.alpha = Q/(Material.eta0*M)+Material.eta0; % Biot coefficient [-]

%% Material parameters (Tian 2023)
% Material.rhof = 952; % fluid density [kg/m3]
% Material.muf = 1e-3; % fluid dynamic viscosity [Pa s]
% Material.Kf = 2.1420e9; % fluid bulk modulus [Pa]
% Material.xif = 2.8e-3; % fluid bulk viscosity [Pa s]
%
% Material.rhos = 2588; % solid density [kg/m3]
% Material.mus = 5.771e9/(1-0.15); % solid shear modulus [Pa]
% Material.Ks = 1.2941e10; % solid bulk modulus [Pa]
%
% Material.rho12 = -404.60; % coupled density [kg/m3]
% Material.alpha = 0.2842; % Biot coefficient [-]
% Material.eta0 = 0.15; % porosity [-]
% Material.k = 1e-13; % permeability [m2] Note: 1D = 1e-12 m2, 1mD = 1e-15 m2

%% Material parameters (Masson, 2006)
% Material.rhof = 1000; % fluid density [kg/m3]
% Material.muf = 1e-3; % fluid dynamic viscosity [Pa s]
% Material.Kf = 2.25e9; % fluid bulk modulus [Pa]
% Material.xif = 0; % fluid bulk viscosity [Pa s]
%
% Material.rhos = 2650; % solid density [kg/m3]
% Material.mus = 455e6/(1-0.3); % solid shear modulus [Pa]
% Material.Ks = 36e9; % solid bulk modulus [Pa]
%
% Material.rho12 = 0; % coupled density [kg/m3]
% Material.alpha = 1-621e6/36e9; % Biot coefficient [-]
% Material.eta0 = 0.3; % porosity [-]
% Material.k = 1e-12; % permeability [m2] Note: 1D = 1e-12 m2, 1mD = 1e-15 m2

%% Material parameters (Ba, 2011)
% Material.rhof = 1040; % fluid density [kg/m3]
% Material.muf = 1e-3; % fluid dynamic viscosity [Pa s]
% Material.Kf = 2.5e9; % fluid bulk modulus [Pa]
% Material.xif = 0; % fluid bulk viscosity [Pa s]
%
% Material.rhos = 2650; % solid density [kg/m3]
% Material.mus = 44e9; % solid shear modulus [Pa]
% Material.Ks = 36e9; % solid bulk modulus [Pa]
%
% Material.rho12 = 0; % coupled density [kg/m3]
% Material.alpha = 1; % Biot coefficient [-]
% Material.eta0 = 0.3; % porosity [-]
% Material.k = 1e-12; % permeability [m2] Note: 1D = 1e-12 m2, 1mD = 1e-15 m2

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

%% Parameter to compare
% permeability [m2]
% mat = logspace(-13, -8);

% porosity
% eta0_min = 0.15;
% eta0_max = 0.45;
% mat = linspace(eta0_min, eta0_max, 20);

% dynamic viscosity [Pa s]
% mat = logspace(-3, 3);

% bulk viscosity [Pa s]
mat = logspace(-3, 3);

%% Initialize variables
vp1_dCS = zeros(length(w), length(mat));
attp1_dCS = zeros(length(w), length(mat));

vp2_dCS = zeros(length(w), length(mat));
attp2_dCS = zeros(length(w), length(mat));

vs1_dCS = zeros(length(w), length(mat));
atts1_dCS = zeros(length(w), length(mat));

vs2_dCS = zeros(length(w), length(mat));
atts2_dCS = zeros(length(w), length(mat));

%% Loop for material parameter range
for i = 1:length(mat)
    % parameter
%     Material.k = mat(i);
%     Material.eta0 = mat(i);
%     Material.muf = mat(i);
    Material.xif = mat(i);
    
    % compute velocity and attenuation
    [vp, attp, vs, atts] = ComputeVelAtt(w, Material);
    % collect data for dCS model
    % fast P
    vp1_dCS(:,i) = vp(:,6);
    attp1_dCS(:,i) = attp(:,6);
    % slow P
    vp2_dCS(:,i) = vp(:,5);
    attp2_dCS(:,i) = attp(:,5);
    % fast S
    vs1_dCS(:,i) = vs(:,6);
    atts1_dCS(:,i) = atts(:,6);
    % slow S
    vs2_dCS(:,i) = vs(:,5);
    atts2_dCS(:,i) = atts(:,5);
end

% contour plots
PlotVelAtt_Contour(w, mat, vp1_dCS, vp2_dCS, vs1_dCS, vs2_dCS, attp1_dCS, attp2_dCS, atts1_dCS, atts2_dCS);
