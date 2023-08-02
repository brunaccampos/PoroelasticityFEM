function PlotConvergenceRate()
% ------------------------------------------------------------------------
% Plot convergence error and slope curve for displacement and pressure
% fields
% ------------------------------------------------------------------------

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
load Results_m1ManBiot_dt1e-5.mat
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
load Results_m2ManBiot_dt1e-5.mat
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
load Results_m3ManBiot_dt1e-5.mat
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
load Results_m4ManBiot_dt1e-5.mat
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
load Results_m5ManBiot_dt1e-5.mat
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
pL2u = polyfit(log(h), log(eL2u),1);
m_L2u = pL2u(1);

pL2p = polyfit(log(h), log(eL2p),1);
m_L2p = pL2p(1);

% Determine slope of energy norm
pENu = polyfit(log(h), log(eENu),1);
m_ENu = pENu(1);

pENp = polyfit(log(h), log(eENp),1);
m_ENp = pENp(1);

% Determine slope of H1 norm
pH1u = polyfit(log(h), log(eH1u),1);
m_H1u = pH1u(1);

pH1p = polyfit(log(h), log(eH1p),1);
m_H1p = pH1p(1);

%% Plots
figure;
% L2 norm displacement
subplot(2,3,1);
loglog(h,eL2u,'g*', 'LineWidth', 1.5);
hold on
loglog(h,exp(pL2u(2))*h.^pL2u(1),'g', 'LineWidth', 1.5);
xlabel('Mesh size (m)');
ylabel('L2-norm');
title(sprintf('L2 - u Convergence: %.2f', m_L2u));
hold off

% L2 norm pressure
subplot(2,3,4);
loglog(h,eL2p,'m*', 'LineWidth', 1.5);
hold on
loglog(h,exp(pL2p(2))*h.^pL2p(1),'m', 'LineWidth', 1.5);
xlabel('Mesh size (m)');
ylabel('L2-norm');
title(sprintf('L2 - p Convergence: %.2f', m_L2p));
hold off

% Energy norm displacement
subplot(2,3,2);
loglog(h,eENu,'b*', 'LineWidth', 1.5);
hold on
loglog(h,exp(pENu(2))*h.^pENu(1),'b', 'LineWidth', 1.5);
xlabel('Mesh size (m)');
ylabel('e-norm');
title(sprintf('Energy - u Convergence: %.2f', m_ENu));
hold off

% Energy norm pressure
subplot(2,3,5);
loglog(h,eENp,'k*', 'LineWidth', 1.5);
hold on
loglog(h,exp(pENp(2))*h.^pENp(1),'k', 'LineWidth', 1.5);
xlabel('Mesh size (m)');
ylabel('q-norm');
title(sprintf('Energy - p Convergence: %.2f', m_ENp));
hold off

% H1 norm displacement
subplot(2,3,3);
loglog(h,eH1u,'r*', 'LineWidth', 1.5);
hold on
loglog(h,exp(pH1u(2))*h.^pH1u(1),'r', 'LineWidth', 1.5);
xlabel('Mesh size (m)');
ylabel('e-norm');
title(sprintf('H1 - u Convergence: %.2f', m_H1u));
hold off

% H1 norm norm pressure
subplot(2,3,6);
loglog(h,eH1p,'c*', 'LineWidth', 1.5);
hold on
loglog(h,exp(pH1p(2))*h.^pH1p(1),'c', 'LineWidth', 1.5);
xlabel('Mesh size (m)');
ylabel('q-norm');
title(sprintf('H1 - p Convergence: %.2f', m_H1p));
hold off 
end