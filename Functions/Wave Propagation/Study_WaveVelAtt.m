% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later
% ------------------------------------------------------------------------
% Study of wave velocity and attenuation: computation for different
% frequency values
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

% uncomment the chosen set of material parameters
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

%% Material parameters (Fellah, 2004) - cancellous bone
% Material.rhof = 1000; % fluid density [kg/m3]
% Material.muf = 1e-3; % fluid dynamic viscosity [Pa s]
% Material.Kf = 2.28e9; % fluid bulk modulus [Pa]
% Material.xif = 3*Material.muf; % fluid bulk viscosity [Pa s]
% 
% Material.rhos = 1960; % solid density [kg/m3]
% Material.mus = 2.6e9; % solid shear modulus [Pa]
% Material.Ks = 20e9; % solid bulk modulus [Pa]
% 
% Material.rho12 = 0; % coupled density [kg/m3]
% Material.alpha = 0.84; % Biot coefficient [-]
% Material.eta0 = 0.83; % porosity [-]
% Material.k = 3e-8; % permeability [m2] Note: 1D = 1e-12 m2, 1mD = 1e-15 m2

%% Material parameters (Hosokawa, 1996) - cancellous bone
% Material.rhof = 930; % fluid density [kg/m3]
% Material.muf = 1.5; % fluid dynamic viscosity [Pa s]
% Material.Kf = 2e9; % fluid bulk modulus [Pa]
% Material.xif = 3*Material.muf; % fluid bulk viscosity [Pa s]
% 
% Material.rhos = 1960; % solid density [kg/m3]
% Material.mus = 2.6e9; % solid shear modulus [Pa]
% Material.Ks = 22e9; % solid bulk modulus [Pa]
% 
% Material.rho12 = 0; % coupled density [kg/m3]
% Material.alpha = 1; % Biot coefficient [-]
% Material.eta0 = 0.8; % porosity [-]
% Material.k = 7e-9; % permeability [m2] Note: 1D = 1e-12 m2, 1mD = 1e-15 m2

%% Material parameters (Cowin, 1999) - cancellous bone
% Material.rhof = 930; % fluid density [kg/m3]
% Material.muf = 1e-3; % fluid dynamic viscosity [Pa s]
% Material.Kf = 2.3e9; % fluid bulk modulus [Pa]
% Material.xif = 3*Material.muf; % fluid bulk viscosity [Pa s]
% 
% Material.rhos = 1960; % solid density [kg/m3]
% Material.mus = 5e9; % solid shear modulus [Pa]
% Material.Ks = 14e9; % solid bulk modulus [Pa]
% 
% Material.rho12 = 0; % coupled density [kg/m3]
% Material.alpha = 0.14; % Biot coefficient [-]
% Material.eta0 = 0.05; % porosity [-]
% Material.k = 1.5e-20; % permeability [m2] Note: 1D = 1e-12 m2, 1mD = 1e-15 m2

%% Material parameters (Picotti, 2024) - Antarctic firn
% Material.rhof = 200; % fluid density [kg/m3]
% Material.muf = 0.1; % fluid dynamic viscosity [Pa s]
% Material.Kf = 571e6; % fluid bulk modulus [Pa]
% Material.xif = Material.muf*3; % fluid bulk viscosity [Pa s]
% 
% Material.rhos = 917; % solid density [kg/m3]
% Material.mus = 5e9; % solid shear modulus [Pa]
% Material.Ks = 10e9; % solid bulk modulus [Pa]
% C = 0.012; 
% Material.rho12 = 0; % coupled density [kg/m3]
% Material.eta0 = 0.80; % porosity [-]
% Material.k = C*Material.eta0^3/Material.rhos^2/(1-Material.eta0)^2; % permeability [m2] Note: 1D = 1e-12 m2, 1mD = 1e-15 m2
% Km = Material.Ks*(1-Material.eta0)^(30.85/(7.76-Material.eta0));
% Material.alpha = 1 - Km/Material.Ks; % Biot coefficient [-]

%% Material parameters (Kimura, 2007) - marine sediments
% Material.rhof = 1025; % fluid density [kg/m3]
% Material.muf = 1e-3; % fluid dynamic viscosity [Pa s]
% Material.Kf = 2.37e9; % fluid bulk modulus [Pa]
% Material.xif = 2.8e-3; % fluid bulk viscosity [Pa s]
% 
% Material.rhos = 2650; % solid density [kg/m3]
% Material.mus = 4.37e5; % solid shear modulus [Pa]
% Material.Ks = 36e9; % solid bulk modulus [Pa]
% 
% Material.rho12 = 0; % coupled density [kg/m3]
% Material.alpha = 1; % Biot coefficient [-]
% Material.eta0 = 0.393; % porosity [-]
% Material.k = 1.15e-10; % permeability [m2] Note: 1D = 1e-12 m2, 1mD = 1e-15 m2

%% Material parameters (Kimura, 2007) - marine sediments and heavy oil
% Material.rhof = 1025; % fluid density [kg/m3]
% Material.muf = 1e1; % fluid dynamic viscosity [Pa s]
% Material.Kf = 2.37e9; % fluid bulk modulus [Pa]
% Material.xif = Material.muf*3; % fluid bulk viscosity [Pa s]
% 
% Material.rhos = 2650; % solid density [kg/m3]
% Material.mus = 4.37e5; % solid shear modulus [Pa]
% Material.Ks = 36e9; % solid bulk modulus [Pa]
% 
% Material.rho12 = 0; % coupled density [kg/m3]
% Material.alpha = 1; % Biot coefficient [-]
% Material.eta0 = 0.393; % porosity [-]
% Material.k = 1.15e-10; % permeability [m2] Note: 1D = 1e-12 m2, 1mD = 1e-15 m2

%% Material parameters (Tang, 2012) - shale
% Material.rhof = 1000; % fluid density [kg/m3]
% Material.muf = 1e-3; % fluid dynamic viscosity [Pa s]
% Material.Kf = 3.3e9; % fluid bulk modulus [Pa]
% Material.xif = 2.8e-3; % fluid bulk viscosity [Pa s]
% 
% Material.rhos = 2600; % solid density [kg/m3]
% Material.mus = 18.4e9; % solid shear modulus [Pa]
% Material.Ks = 35e9; % solid bulk modulus [Pa]
% 
% Material.rho12 = 0; % coupled density [kg/m3]
% Material.alpha = 0.47; % Biot coefficient [-]
% Material.eta0 = 0.08; % porosity [-]
% Material.k = 1e-21; % permeability [m2] Note: 1D = 1e-12 m2, 1mD = 1e-15 m2

%% Material parameters - dCS model
% micro heterogeneity coefficient [-] (Quiroga, 2007)
Material.c = 0;
% Material.c = -0.5;
% Material.c = -0.9;

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

%% Critical frequency limit
wc = Material.muf*Material.eta0/Material.k/Material.rhof;
fprintf('Critical frequency limit is %.2d Hz \n', wc);

%% Compute velocity and attenuation
[vp, attp, vs, atts] = ComputeVelAtt(w, Material);

% plots frequency vs velocity/attenuation
PlotVelAtt(w, vp, attp, vs, atts, wc);

%% Porosity-pressure waves
[vpPorPres, attpPorPres] = ComputeVelAtt_PorPres(w, Material);

% plots frequency vs velocity/attenuation
PlotVelAtt_PorPres(w, vpPorPres, attpPorPres, wc);
