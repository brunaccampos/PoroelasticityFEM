function PlotConvergenceRate()
% ------------------------------------------------------------------------
% Plot convergence error and slope curve for displacement and pressure
% fields
% ------------------------------------------------------------------------

clearvars
clear, clc
close all

%% Initialize variables
nsims = 4;
h = zeros(nsims,1);
eL2u = zeros(nsims,1);
eL2p = zeros(nsims,1);

%% Mesh 1
% load file
load Results_m1_space_Boone.mat
% mesh size
h(1) = ErrorComp.h;
% error displacement
eL2u(1) = ErrorComp.eL2u;
% error pressure
eL2p(1) = ErrorComp.eL2p;

%% Mesh 2
load Results_m2_space_Boone.mat
% mesh size
h(2) = ErrorComp.h;
% error displacement
eL2u(2) = ErrorComp.eL2u;
% error pressure
eL2p(2) = ErrorComp.eL2p;

%% Mesh 3
load Results_m3_space_Boone.mat
% mesh size
h(3) = ErrorComp.h;
% error displacement
eL2u(3) = ErrorComp.eL2u;
% error pressure
eL2p(3) = ErrorComp.eL2p;

%% Mesh 4
load Results_m4_space_Boone.mat
% mesh size
h(4) = ErrorComp.h;
% error displacement
eL2u(4) = ErrorComp.eL2u;
% error pressure
eL2p(4) = ErrorComp.eL2p;

%% Mesh 5
% load Results_m5_space_Boone.mat
% % mesh size
% h(5) = ErrorComp.h;
% % error displacement
% eL2u(5) = ErrorComp.eL2u;
% % error pressure
% eL2p(5) = ErrorComp.eL2p;

%% Curve slope
% Determine slope of L2 norm
pL2u = polyfit(log(h), log(eL2u),1);
m_L2u = pL2u(1);

pL2p = polyfit(log(h), log(eL2p),1);
m_L2p = pL2p(1);

%% Plots
% displacement
figure;
loglog(h,eL2u,'g-o', 'LineWidth', 1.5);
xlabel('Mesh size (m)');
ylabel('L2-norm');
title(sprintf('Convergence for displacement order %.2f', m_L2u));

% pressure
figure;
loglog(h,eL2p,'m-o', 'LineWidth', 1.5);
xlabel('Mesh size (m)');
ylabel('L2-norm');
title(sprintf('Convergence for pressure order %.2f', m_L2p));

end