function PlotVelAtt_Contour(w, mat, vp1_dCS, vp2_dCS, vs1_dCS, vs2_dCS, attp1_dCS, attp2_dCS, atts1_dCS, atts2_dCS)
% Contour plots for frequency vs velocity/attenuation vs material parameter
% for dCS theory
% ------------------------------------------------------------------------

% initialize figure
figure;
tiledlayout(2,4);

%% Plot velocity
% fast P wave
nexttile;
contourf(mat, w, vp1_dCS, 20);
% surf(mat, w, vp1_dCS);
hold on
grid on
title('Fast P wave');
xlabel('Material parameter');
ylabel('Frequency [Hz]');
legend('Velocity [m/s]');
set(gca, 'XScale', 'log');
hold off

% slow P wave
nexttile;
contourf(mat, w, vp2_dCS, 20);
% surf(mat, w, vp2_dCS);
hold on
grid on
title('Slow P wave');
xlabel('Material parameter');
ylabel('Frequency [Hz]');
legend('Velocity [m/s]');
set(gca, 'XScale', 'log');
hold off

% fast S wave
nexttile;
contourf(mat, w, vs1_dCS, 20);
% surf(mat, w, vs1_dCS);
hold on
grid on
title('Fast S wave');
xlabel('Material parameter');
ylabel('Frequency [Hz]');
legend('Velocity [m/s]');
set(gca, 'XScale', 'log');
hold off

% slow S wave
nexttile;
contourf(mat, w, vs2_dCS, 20);
% surf(mat, w, vs2_dCS);
hold on
grid on
title('Slow S wave');
xlabel('Material parameter');
ylabel('Frequency [Hz]');
legend('Velocity [m/s]');
set(gca, 'XScale', 'log');
hold off

%% Plot attenuation
% fast P wave
nexttile;
contourf(mat, w, attp1_dCS, 20);
% surf(mat, w, attp1_dCS);
hold on
grid on
title('Fast P wave');
xlabel('Material parameter');
ylabel('Frequency [Hz]');
legend('Attenuation [1/m]');
set(gca, 'XScale', 'log');
hold off

% slow P wave
nexttile;
contourf(mat, w, attp2_dCS, 20);
% surf(mat, w, attp2_dCS);
hold on
grid on
title('Slow P wave');
xlabel('Material parameter');
ylabel('Frequency [Hz]');
legend('Attenuation [1/m]');
set(gca, 'XScale', 'log');
hold off

% fast S wave
nexttile;
contourf(mat, w, atts1_dCS, 20);
% surf(mat, w, atts1_dCS);
hold on
grid on
title('Fast S wave');
xlabel('Material parameter');
ylabel('Frequency [Hz]');
legend('Attenuation [1/m]');
set(gca, 'XScale', 'log');
hold off

% slow S wave
nexttile;
contourf(mat, w, atts2_dCS, 20);
% surf(mat, w, atts2_dCS);
hold on
grid on
title('Slow S wave');
xlabel('Material parameter');
ylabel('Frequency [Hz]');
legend('Attenuation [1/m]');
set(gca, 'XScale', 'log');
hold off

end