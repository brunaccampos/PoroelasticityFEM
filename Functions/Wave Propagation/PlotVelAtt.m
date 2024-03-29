function PlotVelAtt(w, vp, attp, vs, atts)
% Plot velocity and attenuation for compressional and shear waves using
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

plotZH = 0;

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
semilogx(w,vp(:,6),'r:', 'LineWidth', 1.5); % dCS
legend('BT', 'dCS');
if plotZH
    semilogx(w,vp(:,4),'b--', 'LineWidth', 1.5); % ZH
    legend('BT', 'dCS', 'ZH');
end
hold off

% second P wave
nexttile;
semilogx(w,vp(:,1), 'k-', 'LineWidth', 1.5); % BT
hold on
grid on
xlabel('Frequency [Hz]');
ylabel('Velocity [m/s]');
title('Second P wave');
semilogx(w,vp(:,5), 'r:', 'LineWidth', 1.5); % dCS
legend('BT', 'dCS');
if plotZH
    semilogx(w,vp(:,3), 'b--', 'LineWidth', 1.5); % ZH
    legend('BT', 'dCS', 'ZH');
end
hold off

% first S wave
nexttile;
semilogx(w,vs(:,1), 'k-', 'LineWidth', 1.5); % BT
hold on
grid on
xlabel('Frequency [Hz]');
ylabel('Velocity [m/s]');
title('First S wave');
semilogx(w,vs(:,6), 'r:', 'LineWidth', 1.5); % dCS
legend('BT', 'dCS');
if plotZH
    semilogx(w,vs(:,4), 'b--', 'LineWidth', 1.5); % ZH
    legend('BT', 'dCS', 'ZH');
end
hold off

% second S wave
nexttile;
semilogx(w,vs(:,5), 'r:', 'LineWidth', 1.5); % dCS
hold on
grid on
xlabel('Frequency [Hz]');
ylabel('Velocity [m/s]');
title('Second S wave');
legend('dCS');
if plotZH
    semilogx(w,vs(:,3), 'b--', 'LineWidth', 1.5); % ZH
    legend('dCS', 'ZH');
end
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
semilogx(w,attp(:,6),'r:', 'LineWidth', 1.5); % dCS
legend('BT', 'dCS');
if plotZH
    semilogx(w,attp(:,4),'b--', 'LineWidth', 1.5); % ZH
    legend('BT', 'dCS', 'ZH');
end
hold off

% second P wave
nexttile;
semilogx(w,attp(:,1), 'k-', 'LineWidth', 1.5); % BT
hold on
grid on
xlabel('Frequency [Hz]');
ylabel('Attenuation');
title('Second P wave');
semilogx(w,attp(:,5), 'r:', 'LineWidth', 1.5); % dCS
legend('BT', 'dCS');
if plotZH
    semilogx(w,attp(:,3), 'b--', 'LineWidth', 1.5); % ZH
    legend('BT', 'dCS', 'ZH');
end
hold off

% first S wave
nexttile;
semilogx(w,atts(:,1), 'k-', 'LineWidth', 1.5); % BT
hold on
grid on
xlabel('Frequency [Hz]');
ylabel('Attenuation');
title('First S wave');
semilogx(w,atts(:,6), 'r:', 'LineWidth', 1.5); % dCS
legend('BT', 'dCS');
if plotZH
    semilogx(w,atts(:,4), 'b--', 'LineWidth', 1.5); % ZH
    legend('BT', 'dCS', 'ZH');
end
hold off

% second S wave
nexttile;
semilogx(w,atts(:,5), 'r:', 'LineWidth', 1.5); % dCS
hold on
grid on
xlabel('Frequency [Hz]');
ylabel('Attenuation');
title('Second S wave');
legend('dCS');
if plotZH
    semilogx(w,atts(:,3), 'b--', 'LineWidth', 1.5); % ZH
    legend('ZH', 'dCS');
end
hold off


end