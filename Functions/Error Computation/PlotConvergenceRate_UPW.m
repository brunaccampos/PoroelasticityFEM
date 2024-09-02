function PlotConvergenceRate_UPW()
% ------------------------------------------------------------------------
% Plot convergence error and slope curve for displacement and pressure
% fields
% ------------------------------------------------------------------------

%% Initialize variables
% number of simulations
nsims = 5;
% mesh size vector
h = zeros(nsims,1);
% L2 norm
eL2us = zeros(nsims,1);
eL2p = zeros(nsims,1);
eL2w = zeros(nsims,1);

%% Mesh 1
% load file
load Results_m11ManUPUBiot.mat
% mesh size
h(1) = ErrorComp.h;
% L2 solid displacement
eL2us(1) = ErrorComp.eL2u;
% L2 pressure
eL2p(1) = ErrorComp.eL2p;
% L2 fluid displacement
eL2w(1) = ErrorComp.eL2w;

%% Mesh 2
load Results_m12ManUPUBiot.mat
% mesh size
h(2) = ErrorComp.h;
% L2 solid displacement
eL2us(2) = ErrorComp.eL2u;
% L2 pressure
eL2p(2) = ErrorComp.eL2p;
% L2 fluid displacement
eL2w(2) = ErrorComp.eL2w;

%% Mesh 3
load Results_m13ManUPUBiot.mat
% mesh size
h(3) = ErrorComp.h;
% L2 solid displacement
eL2us(3) = ErrorComp.eL2u;
% L2 pressure
eL2p(3) = ErrorComp.eL2p;
% L2 fluid displacement
eL2w(3) = ErrorComp.eL2w;

%% Mesh 4
load Results_m14ManUPUBiot.mat
% mesh size
h(4) = ErrorComp.h;
% L2 solid displacement
eL2us(4) = ErrorComp.eL2u;
% L2 pressure
eL2p(4) = ErrorComp.eL2p;
% L2 fluid displacement
eL2w(4) = ErrorComp.eL2w;

%% Mesh 5
load Results_m15ManUPUBiot.mat
% mesh size
h(5) = ErrorComp.h;
% L2 solid displacement
eL2us(5) = ErrorComp.eL2u;
% L2 pressure
eL2p(5) = ErrorComp.eL2p;
% L2 fluid displacement
eL2w(5) = ErrorComp.eL2w;

%% Curve slope
% Determine slope of L2 norm
pL2us = polyfit(log(h), log(eL2us),1);
m_L2us = pL2us(1);

pL2p = polyfit(log(h), log(eL2p),1);
m_L2p = pL2p(1);

pL2w = polyfit(log(h), log(eL2w),1);
m_L2w = pL2w(1);

%% Plots
figure;
tiledlayout(1,3);
% L2 norm solid displacement
nexttile;
loglog(h,eL2us,'g*', 'LineWidth', 1.5);
hold on
loglog(h,exp(pL2us(2))*h.^pL2us(1),'g', 'LineWidth', 1.5);
grid on
xlabel('Mesh size (m)');
ylabel('L2-norm');
title(sprintf('L2 - us Convergence: %.2f', m_L2us));
hold off

% L2 norm pressure
nexttile;
loglog(h,eL2p,'m*', 'LineWidth', 1.5);
hold on
loglog(h,exp(pL2p(2))*h.^pL2p(1),'m', 'LineWidth', 1.5);
grid on
xlabel('Mesh size (m)');
ylabel('L2-norm');
title(sprintf('L2 - p Convergence: %.2f', m_L2p));
hold off

% L2 norm relative fluid velocity
nexttile;
loglog(h,eL2w,'r*', 'LineWidth', 1.5);
hold on
loglog(h,exp(pL2w(2))*h.^pL2w(1),'r', 'LineWidth', 1.5);
grid on
xlabel('Mesh size (m)');
ylabel('L2-norm');
title(sprintf('L2 - w Convergence: %.2f', m_L2w));
hold off

end