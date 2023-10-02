function PlotNDPressureEqCoef(Material)
% ------------------------------------------------------------------------
% Plot non dimensional coefficients from pressure equation to study the
% effect of the parameter n
% ------------------------------------------------------------------------

% parameter n
n = (0:0.01:Material.Ks/Material.Kf);

% deltas
deltaF = (Material.alpha - Material.n) * Material.n * n / Material.Ks ./ (Material.Minv - (1-n)*(Material.alpha - Material.n)/Material.Ks);
deltaS = (Material.alpha - Material.n) * Material.n / Material.Kf ./ (Material.Minv - (1-n)*(Material.alpha - Material.n)/Material.Ks);

% find position where n = 1
index = find (n == 1);

% coefficient for Kps
hi1 = deltaS * Material.Kf / (Material.n * Material.E);
% values normalized with n=1
hi1norm = hi1./hi1(index);

% coefficient for Kpf
hi2 = (1 - deltaF/Material.n) * Material.Kf / Material.E;
% normalized values
hi2norm = hi2./hi2(index);

figure;
% plot regular values
subplot(1,2,1);
plot(n, hi1, 'LineWidth', 1.5);
hold on
plot(n, hi2, 'LineWidth', 1.5);
grid on
xline(1, '--', 'LineWidth', 1.5);
xlabel('n');
ylabel('\zeta');
title('Variation of \zeta_1 and \zeta_2');
legend('\zeta_1', '\zeta_2', 'n=1');
hold off

% plot normalized values
subplot(1,2,2);
plot(n, hi1norm, 'LineWidth', 1.5);
hold on
plot(n, hi2norm, 'LineWidth', 1.5);
grid on
xline(1, '--', 'LineWidth', 1.5);
xlabel('n');
ylabel('\zeta / \zeta(n=1)');
title('Variation of \zeta_1 and \zeta_2');
legend('\zeta_1', '\zeta_2', 'n=1');
hold off

end