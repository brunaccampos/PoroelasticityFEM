function PlotVelAtt_PorPres(w, vp, attp, wc)
% Plot velocity and attenuation for porosity-pressure waves
% ------------------------------------------------------------------------

% initialize figure
figure;
tiledlayout(2,2);

%% Plots phase velocity
% wave 1
nexttile;
semilogx(w,vp(:,1),'b-', 'LineWidth', 1.5);
hold on
grid on
xline(wc, '--k', 'LineWidth', 1.5); % critical frequency
xlabel('Frequency [Hz]');
ylabel('Velocity [m/s]');
title('PorPre - Slow P');
hold off

% wave 3
nexttile;
semilogx(w,vp(:,3),'b-', 'LineWidth', 1.5);
hold on
grid on
xline(wc, '--k', 'LineWidth', 1.5); % critical frequency
xlabel('Frequency [Hz]');
ylabel('Velocity [m/s]');
title('PorPre - Fast P');
hold off

%% Plots attenuation
% wave 1
nexttile;
semilogx(w,attp(:,1),'m-', 'LineWidth', 1.5);
hold on
grid on
xline(wc, '--k', 'LineWidth', 1.5); % critical frequency
xlabel('Frequency [Hz]');
ylabel('Attenuation');
title('PorPre - Slow P');
hold off

% wave 3
nexttile;
semilogx(w,attp(:,3),'m-', 'LineWidth', 1.5);
hold on
grid on
xline(wc, '--k', 'LineWidth', 1.5); % critical frequency
xlabel('Frequency [Hz]');
ylabel('Attenuation');
title('PorPre - Fast P');
hold off

end