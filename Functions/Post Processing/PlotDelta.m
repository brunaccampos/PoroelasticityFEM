function PlotDelta(Material)
% ------------------------------------------------------------------------
% Plot variation of deltaF and deltaS for different values of n
% ------------------------------------------------------------------------

% parameter n
n = (0:0.01:Material.Ks/Material.Kf);
% deltas
deltaF = (Material.alpha - Material.n) * Material.n * n / Material.Ks ./ (Material.Minv - (1-n)*(Material.alpha - Material.n)/Material.Ks);
deltaS = (Material.alpha - Material.n) * Material.n / Material.Kf ./ (Material.Minv - (1-n)*(Material.alpha - Material.n)/Material.Ks);

% plot
figure;
plot(n, deltaF, 'LineWidth', 1.5);
hold on
plot(n, deltaS, 'LineWidth', 1.5);
xline(1, '--', 'LineWidth', 1.5);
xlabel('n');
ylabel('\delta');
title('Variation of \delta_s and \delta_f');
legend('\delta_f', '\delta_s', 'n=1');

end