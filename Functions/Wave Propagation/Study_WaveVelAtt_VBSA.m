% Study of wave velocity and attenuation: variance-based sensitivity
% analysis
% ------------------------------------------------------------------------

clear vars
clc
close all
fprintf('Variance-based sensitivity analysis: wave velocity and attenuation \n');
fprintf('Define parameters ...\n');

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

%% Wave mode to compute
% type = 'Fast P - Velocity'; index = 6; location1 = [0.7648 0.5 0.1210 0.1719]; location2 = [0.62 0.51 0.1 0.1];
% type = 'Slow P - Velocity'; index = 5;  location1 = [0.7648 0.5 0.1210 0.1719]; location2 = [0.62 0.51 0.1 0.1];
% type = 'Fast S - Velocity'; index = 6;  location1 = [0.7648 0.5 0.1210 0.1719]; location2 = [0.62 0.51 0.1 0.1];
% type = 'Slow S - Velocity'; index = 5;  location1 = [0.7648 0.5 0.1210 0.1719]; location2 = [0.62 0.51 0.1 0.1];
% type = 'Fast P - Attenuation'; index = 6; location1 = [0.7648 0.5 0.1210 0.1719]; location2 = [0.62 0.51 0.1 0.1];
% type = 'Slow P - Attenuation'; index = 5; location1 = [0.7648 0.45 0.1210 0.1719]; location2 = [0.62 0.46 0.1 0.1];
% type = 'Fast S - Attenuation'; index = 6; location1 = [0.7648 0.5 0.1210 0.1719]; location2 = [0.62 0.51 0.1 0.1];
type = 'Slow S - Attenuation'; index = 5; location1 = [0.7648 0.5 0.1210 0.1719]; location2 = [0.62 0.51 0.1 0.1];

%% Frequency array
w = linspace(1e0, 1e6, 10); % frequency [Hz]

%% Parameters to compare
% permeability [m2]
mat_k = linspace(1e-21, 1e-13, 10);

% porosity
eta0_min = 0.05;
eta0_max = 0.35;
mat_eta0 = linspace(eta0_min, eta0_max, 10);

% dynamic viscosity/ bulk viscosity [Pa s]
mat_muf = linspace(1e-3, 1e-1, 10);

% dimensionless parameter
mat_c = linspace(-1, 0, 5);

% porosity effective pressure coefficient
mat_n = linspace(0, Material.Ks/Material.Kf, 10);

%% Sobol analysis parameters
S_first = zeros(length(w), 5);
S_total = zeros(length(w), 5);
Y = zeros(length(w), length(mat_eta0), length(mat_k), length(mat_muf), length(mat_c), length(mat_n));

%% Loop for material parameter range
fprintf('Loop over material parameter ranges...\n');
for ieta = 1:length(mat_eta0)
    for ik = 1:length(mat_k)
        for imu = 1:length(mat_muf)
            for ic = 1:length(mat_c)
                for in = 1:length(mat_n)
                    % parameters
                    Material.eta0 = mat_eta0(ieta);
                    Material.k = mat_k(ik);
                    Material.muf = mat_muf(imu);
                    Material.c = mat_c(ic);
                    Material.n = mat_n(in);

                    % compute velocity and attenuation
                    [vp, attp, vs, atts] = ComputeVelAtt(w, Material);

                    % pick response and mode
                    switch type
                        case 'Fast P - Velocity'
                            Y(:, ieta, ik, imu, ic, in) = vp(:,index);
                        case 'Fast P - Attenuation'
                            Y(:, ieta, ik, imu, ic, in) = attp(:,index);
                        case 'Fast S - Velocity'
                            Y(:, ieta, ik, imu, ic, in) = vs(:,index);
                        case 'Fast S - Attenuation'
                            Y(:, ieta, ik, imu, ic, in) = atts(:,index);
                        case 'Slow P - Velocity'
                            Y(:, ieta, ik, imu, ic, in) = vp(:,index);
                        case 'Slow P - Attenuation'
                            Y(:, ieta, ik, imu, ic, in) = attp(:,index);
                        case 'Slow S - Velocity'
                            Y(:, ieta, ik, imu, ic, in) = vs(:,index);
                        case 'Slow S - Attenuation'
                            Y(:, ieta, ik, imu, ic, in) = atts(:,index);
                    end
                end
            end
        end
    end
end

%% Loop over frequencies
fprintf('Compute Sobol indices...\n');
for i = 1:length(w)
    Yf = squeeze(Y(i,:,:,:,:,:));
    % total variance (flattened)
    VarY = var(Yf(:), 1);  % use population variance (normalization by N)
    if VarY == 0
        S_first(i,:) = 0;
        S_total(i,:) = 0;
        continue;
    end

    % First order indices: Var(E[Y | Xi])/Var(Y)
    % E[Y | eta0]: average over k, mu, c, n
    E_Y_given_eta0 = squeeze(mean(mean(mean(mean(Yf, 5), 4), 3), 2));
    V_eta0 = var(E_Y_given_eta0, 1);

    % E[Y | k]: average over eta0, mu, c, n
    E_Y_given_k = squeeze(mean(mean(mean(mean(Yf, 5), 4), 3), 1));
    V_k = var(E_Y_given_k, 1);

    % E[Y | mu]: average over eta0, k, c, n
    E_Y_given_muf = squeeze(mean(mean(mean(mean(Yf, 5), 4), 2), 1));
    V_muf = var(E_Y_given_muf, 1);

    % E[Y | c]: average over eta0, k, mu, n
    E_Y_given_c = squeeze(mean(mean(mean(mean(Yf, 5), 3), 2), 1));
    V_c = var(E_Y_given_c, 1);

    % E[Y | n]: average over eta0, k, mu, c
    E_Y_given_n = squeeze(mean(mean(mean(mean(Yf, 4), 3), 2), 1));
    V_n = var(E_Y_given_n, 1);

    S_first(i, :) = [V_eta0, V_k, V_muf, V_c, V_n] / VarY;

    % Total order indices: ST_i = 1 - Var(E[Y | X_~i])/Var(Y)
    % E[Y | not eta0]: average over eta0
    E_Y_given_noteta0 = squeeze(mean(Yf, 1));
    V_noteta0 = var(E_Y_given_noteta0(:), 1);
    ST_eta0 = 1 - V_noteta0 / VarY;

    % E[Y | not k]: average over k
    E_Y_given_notk = squeeze(mean(Yf, 2));
    V_notk = var(E_Y_given_notk(:), 1);
    ST_k = 1 - V_notk / VarY;

    % E[Y | not mu]: average over mu
    E_Y_given_notmuf = squeeze(mean(Yf, 3));
    V_notmuf = var(E_Y_given_notmuf(:), 1);
    ST_muf = 1 - V_notmuf / VarY;

    % E[Y | not c]: average over c
    E_Y_given_notc = squeeze(mean(Yf, 4));
    V_notc = var(E_Y_given_notc(:), 1);
    ST_c = 1 - V_notc / VarY;

    % E[Y | not n]: average over n
    E_Y_given_notn = squeeze(mean(Yf, 5));
    V_notn = var(E_Y_given_notn(:), 1);
    ST_n = 1 - V_notn / VarY;

    S_total(i, :) = [ST_eta0, ST_k, ST_muf, ST_c, ST_n];
end

fprintf('Generate plots...\n');

% color options
colors = orderedcolors("gem");

% % Number of markers desired
% numMarkers = 30;
% 
% % Create logarithmically spaced x-coordinates for markers
% xmarkers = logspace(log10(w(1)), log10(w(end)), numMarkers);
% 
% % Interpolate y-values at these x-marker locations (in log space for accuracy)
% ymarkers = 10.^interp1(log10(w), log10(S_first(:,1)), log10(xmarkers));


% plot first-order and total-order vs frequency
figure;
hold on
grid minor
xlim tight
ylim padded
title(sprintf('%s', type),'interpreter','latex', 'FontSize', 14);
plot(w, S_first(:,1), '--o', 'Color', colors(1,:));
plot(w, S_first(:,2), '--square', 'Color', colors(2,:));
plot(w, S_first(:,3), '--^', 'Color', colors(3,:));
plot(w, S_first(:,4), '--diamond', 'Color', colors(4,:));
plot(w, S_first(:,5), '--+', 'Color', colors(5,:));
plt(1) = plot(w, S_total(:,1), '-o', 'Color', colors(1,:));
plt(2) = plot(w, S_total(:,2), '-square', 'Color', colors(2,:));
plt(3) = plot(w, S_total(:,3), '-^', 'Color', colors(3,:));
plt(4) = plot(w, S_total(:,4), '-diamond', 'Color', colors(4,:));
plt(5) = plot(w, S_total(:,5), '-+', 'Color', colors(5,:));

% dummy plot for indices legend
plt_first = plot(nan, nan, '--k');
plt_total = plot(nan, nan, '-k');

set(gca,'TickLabelInterpreter','latex', 'FontSize', 14);
xlabel('Frequency [Hz]','interpreter','latex', 'FontSize', 14);
ylabel('Sobol Index','interpreter','latex', 'FontSize', 14);
ylim([0 1]);
leg = legend(plt,'$\eta_0$','$k$','$\mu_f$', '$\delta_{\mu}$', '$n$', 'Location', location1, 'interpreter','latex', 'FontSize', 14);
set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.5);

ax = axes('position',get(gca,'position'),'visible','off');
legend(ax, [plt_first plt_total], 'S1', 'ST', ...
    'Interpreter','latex', ...
    'Location', location2,'FontSize',14);

% save figure
name_png = sprintf('VBSA_%s.png', type);
name_pdf = sprintf('VBSA_%s.pdf', type);
exportgraphics(gcf, name_png, 'Resolution', 600);
exportgraphics(gcf, name_pdf);

fprintf('Plotting and saving images complete!\n');
