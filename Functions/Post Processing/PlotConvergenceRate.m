function PlotConvergenceRate()
% ------------------------------------------------------------------------
% Plot convergence error and slope curve for displacement and pressure
% fields
% ------------------------------------------------------------------------

clearvars
clear, clc
close all

%% Initialize variables
nsims = 5;
h = zeros(nsims,1);
eL2u = zeros(nsims,1);
eL2p = zeros(nsims,1);
eENu = zeros(nsims,1);
eENp = zeros(nsims,1);
eH1u = zeros(nsims,1);
eH1p = zeros(nsims,1);

%% Mesh 1
% load file
load Results_m1Boone_testH1.mat
% mesh size
h(1) = ErrorComp.h;
% error displacement
eL2u(1) = ErrorComp.eL2u;
% error pressure
eL2p(1) = ErrorComp.eL2p;
% error strain
eENu(1) = ErrorComp.eENu;
% error flux
eENp(1) = ErrorComp.eENp;
% error H1
eH1u(1) = ErrorComp.eH1u;
% error H1
eH1p(1) = ErrorComp.eH1p;

%% Mesh 2
load Results_m2Boone_testH1.mat
% mesh size
h(2) = ErrorComp.h;
% error displacement
eL2u(2) = ErrorComp.eL2u;
% error pressure
eL2p(2) = ErrorComp.eL2p;
% error strain
eENu(2) = ErrorComp.eENu;
% error flux
eENp(2) = ErrorComp.eENp;
% error H1
eH1u(2) = ErrorComp.eH1u;
% error H1
eH1p(2) = ErrorComp.eH1p;

%% Mesh 3
load Results_m3Boone_testH1.mat
% mesh size
h(3) = ErrorComp.h;
% error displacement
eL2u(3) = ErrorComp.eL2u;
% error pressure
eL2p(3) = ErrorComp.eL2p;
% error strain
eENu(3) = ErrorComp.eENu;
% error flux
eENp(3) = ErrorComp.eENp;
% error H1
eH1u(3) = ErrorComp.eH1u;
% error H1
eH1p(3) = ErrorComp.eH1p;

%% Mesh 4
load Results_m4Boone_testH1.mat
% mesh size
h(4) = ErrorComp.h;
% error displacement
eL2u(4) = ErrorComp.eL2u;
% error pressure
eL2p(4) = ErrorComp.eL2p;
% error strain
eENu(4) = ErrorComp.eENu;
% error flux
eENp(4) = ErrorComp.eENp;
% error H1
eH1u(4) = ErrorComp.eH1u;
% error H1
eH1p(4) = ErrorComp.eH1p;

%% Mesh 5
load Results_m5Boone_testH1.mat
% mesh size
h(5) = ErrorComp.h;
% error displacement
eL2u(5) = ErrorComp.eL2u;
% error pressure
eL2p(5) = ErrorComp.eL2p;
% error strain
eENu(5) = ErrorComp.eENu;
% error flux
eENp(5) = ErrorComp.eENp;
% error H1
eH1u(5) = ErrorComp.eH1u;
% error H1
eH1p(5) = ErrorComp.eH1p;

%% Curve slope
% Determine slope of L2 norm
pL2u = polyfit(log(h(1:2)), log(eL2u(1:2)),1);
m_L2u = pL2u(1);

pL2p = polyfit(log(h(1:2)), log(eL2p(1:2)),1);
m_L2p = pL2p(1);

% Determine slope of energy norm
pENu = polyfit(log(h(1:2)), log(eENu(1:2)),1);
m_ENu = pENu(1);

pENp = polyfit(log(h(1:2)), log(eENp(1:2)),1);
m_ENp = pENp(1);

% Determine slope of H1 norm
pH1u = polyfit(log(h(1:2)), log(eH1u(1:2)),1);
m_H1u = pH1u(1);

pH1p = polyfit(log(h(1:2)), log(eH1p(1:2)),1);
m_H1p = pH1p(1);

%% Plots
% L2 norm displacement
figure;
loglog(h,eL2u,'g-o', 'LineWidth', 1.5);
xlabel('Mesh size (m)');
ylabel('L2-norm');
title(sprintf('L2 - u Convergence: %.2f', m_L2u));

% L2 norm pressure
figure;
loglog(h,eL2p,'m-o', 'LineWidth', 1.5);
xlabel('Mesh size (m)');
ylabel('L2-norm');
title(sprintf('L2 - p Convergence: %.2f', m_L2p));

% Energy norm displacement
figure;
loglog(h,eENu,'b-o', 'LineWidth', 1.5);
xlabel('Mesh size (m)');
ylabel('e-norm');
title(sprintf('Energy - u Convergence: %.2f', m_ENu));

% Energy norm pressure
figure;
loglog(h,eENp,'k-o', 'LineWidth', 1.5);
xlabel('Mesh size (m)');
ylabel('q-norm');
title(sprintf('Energy - p Convergence: %.2f', m_ENp));

% H1 norm displacement
figure;
loglog(h,eH1u,'r-o', 'LineWidth', 1.5);
xlabel('Mesh size (m)');
ylabel('e-norm');
title(sprintf('H1 - u Convergence: %.2f', m_H1u));

% H1 norm norm pressure
figure;
loglog(h,eH1p,'c-o', 'LineWidth', 1.5);
xlabel('Mesh size (m)');
ylabel('q-norm');
title(sprintf('H1 - p Convergence: %.2f', m_H1p));

end