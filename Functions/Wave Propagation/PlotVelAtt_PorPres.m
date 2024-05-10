function PlotVelAtt_PorPres(w, vp, attp)
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
% initialize figure
figure;
tiledlayout(2,4);

%% Plots phase velocity
% wave 1
nexttile;
semilogx(w,vp(:,1),'b-', 'LineWidth', 1.5);
hold on
grid on
xlabel('Frequency [Hz]');
ylabel('Velocity [m/s]');
title('Wave 1');
hold off

% wave 2
nexttile;
semilogx(w,vp(:,2), 'b-', 'LineWidth', 1.5);
hold on
grid on
xlabel('Frequency [Hz]');
ylabel('Velocity [m/s]');
title('Wave 2');
hold off

% wave 3
nexttile;
semilogx(w,vp(:,3),'b-', 'LineWidth', 1.5);
hold on
grid on
xlabel('Frequency [Hz]');
ylabel('Velocity [m/s]');
title('Wave 3');
hold off

% wave 4
nexttile;
semilogx(w,vp(:,4), 'b-', 'LineWidth', 1.5);
hold on
grid on
xlabel('Frequency [Hz]');
ylabel('Velocity [m/s]');
title('Wave 4');
hold off

%% Plots attenuation
% wave 1
nexttile;
semilogx(w,attp(:,1),'r-', 'LineWidth', 1.5);
hold on
grid on
xlabel('Frequency [Hz]');
ylabel('Attenuation');
title('Wave 1');
hold off

% wave 2
nexttile;
semilogx(w,attp(:,2), 'r-', 'LineWidth', 1.5);
hold on
grid on
xlabel('Frequency [Hz]');
ylabel('Attenuation');
title('Wave 2');
hold off

% wave 3
nexttile;
semilogx(w,attp(:,3),'r-', 'LineWidth', 1.5);
hold on
grid on
xlabel('Frequency [Hz]');
ylabel('Attenuation');
title('Wave 3');
hold off

% wave 4
nexttile;
semilogx(w,attp(:,4), 'r-', 'LineWidth', 1.5);
hold on
grid on
xlabel('Frequency [Hz]');
ylabel('Attenuation');
title('Wave 4');
hold off

end