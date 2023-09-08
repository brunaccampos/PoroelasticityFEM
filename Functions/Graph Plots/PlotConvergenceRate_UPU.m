function PlotConvergenceRate_UPU()
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
eL2uf = zeros(nsims,1);
% Energy norm
eENus = zeros(nsims,1);
eENp = zeros(nsims,1);
eENuf = zeros(nsims,1);

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
eL2uf(1) = ErrorComp.eL2uf;
% EN solid displacement
eENus(1) = ErrorComp.eENu;
% EN pressure
eENp(1) = ErrorComp.eENp;
% EN solid displacement
eENuf(1) = ErrorComp.eENuf;

%% Mesh 2
load Results_m12ManUPUBiot.mat
% mesh size
h(2) = ErrorComp.h;
% L2 solid displacement
eL2us(2) = ErrorComp.eL2u;
% L2 pressure
eL2p(2) = ErrorComp.eL2p;
% L2 fluid displacement
eL2uf(2) = ErrorComp.eL2uf;
% EN solid displacement
eENus(2) = ErrorComp.eENu;
% EN pressure
eENp(2) = ErrorComp.eENp;
% EN solid displacement
eENuf(2) = ErrorComp.eENuf;

%% Mesh 3
load Results_m13ManUPUBiot.mat
% mesh size
h(3) = ErrorComp.h;
% L2 solid displacement
eL2us(3) = ErrorComp.eL2u;
% L2 pressure
eL2p(3) = ErrorComp.eL2p;
% L2 fluid displacement
eL2uf(3) = ErrorComp.eL2uf;
% EN solid displacement
eENus(3) = ErrorComp.eENu;
% EN pressure
eENp(3) = ErrorComp.eENp;
% EN solid displacement
eENuf(3) = ErrorComp.eENuf;

%% Mesh 4
load Results_m14ManUPUBiot.mat
% mesh size
h(4) = ErrorComp.h;
% L2 solid displacement
eL2us(4) = ErrorComp.eL2u;
% L2 pressure
eL2p(4) = ErrorComp.eL2p;
% L2 fluid displacement
eL2uf(4) = ErrorComp.eL2uf;
% EN solid displacement
eENus(4) = ErrorComp.eENu;
% EN pressure
eENp(4) = ErrorComp.eENp;
% EN solid displacement
eENuf(4) = ErrorComp.eENuf;

%% Mesh 5
load Results_m15ManUPUBiot.mat
% mesh size
h(5) = ErrorComp.h;
% L2 solid displacement
eL2us(5) = ErrorComp.eL2u;
% L2 pressure
eL2p(5) = ErrorComp.eL2p;
% L2 fluid displacement
eL2uf(5) = ErrorComp.eL2uf;
% EN solid displacement
eENus(5) = ErrorComp.eENu;
% EN pressure
eENp(5) = ErrorComp.eENp;
% EN solid displacement
eENuf(5) = ErrorComp.eENuf;

%% Curve slope
% Determine slope of L2 norm
pL2us = polyfit(log(h), log(eL2us),1);
m_L2us = pL2us(1);

pL2p = polyfit(log(h), log(eL2p),1);
m_L2p = pL2p(1);

pL2uf = polyfit(log(h), log(eL2uf),1);
m_L2uf = pL2uf(1);

% Determine slope of energy norm
pENus = polyfit(log(h), log(eENus),1);
m_ENus = pENus(1);

pENp = polyfit(log(h), log(eENp),1);
m_ENp = pENp(1);

pENuf = polyfit(log(h), log(eENuf),1);
m_ENuf = pENuf(1);

%% Plots
tiledlayout(2,3);
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

% L2 norm fluid displacement
nexttile;
loglog(h,eL2uf,'r*', 'LineWidth', 1.5);
hold on
loglog(h,exp(pL2uf(2))*h.^pL2uf(1),'r', 'LineWidth', 1.5);
grid on
xlabel('Mesh size (m)');
ylabel('L2-norm');
title(sprintf('L2 - uf Convergence: %.2f', m_L2uf));
hold off

% Energy norm solid displacement
nexttile;
loglog(h,eENus,'b*', 'LineWidth', 1.5);
hold on
loglog(h,exp(pENus(2))*h.^pENus(1),'b', 'LineWidth', 1.5);
grid on
xlabel('Mesh size (m)');
ylabel('EN-norm');
title(sprintf('EN - us Convergence: %.2f', m_ENus));
hold off

% Energy norm pressure
nexttile;
loglog(h,eENp,'k*', 'LineWidth', 1.5);
hold on
loglog(h,exp(pENp(2))*h.^pENp(1),'k', 'LineWidth', 1.5);
grid on
xlabel('Mesh size (m)');
ylabel('EN-norm');
title(sprintf('EN - p Convergence: %.2f', m_ENp));
hold off

% Energy norm fluid displacement
nexttile;
loglog(h,eENuf,'c*', 'LineWidth', 1.5);
hold on
loglog(h,exp(pENuf(2))*h.^pENuf(1),'c', 'LineWidth', 1.5);
grid on
xlabel('Mesh size (m)');
ylabel('EN-norm');
title(sprintf('EN - uf Convergence: %.2f', m_ENuf));
hold off 
end