% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function PlotDelta(Material)
% Plot variation of deltaF and deltaS for different values of n

% parameter n
n = (0:0.01:Material.Ks/Material.Kf);
% deltas
deltaf = (Material.alpha - Material.eta0) * Material.eta0 * n / Material.Ks ./ (Material.Minv - (1-n)*(Material.alpha - Material.eta0)/Material.Ks);
deltas = (Material.alpha - Material.eta0) * Material.eta0 / Material.Kf ./ (Material.Minv - (1-n)*(Material.alpha - Material.eta0)/Material.Ks);

% plot
figure;
plot(n, deltaf, 'LineWidth', 1.5);
hold on
plot(n, deltas, 'LineWidth', 1.5);
grid on
xline(1, '--', 'LineWidth', 1.5);
xlabel('n');
ylabel('\delta');
title('Variation of \delta_s and \delta_f');
legend('\delta_f', '\delta_s', 'n=1');

end