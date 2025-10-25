% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function PlotGraphsUP_OverTime(MeshU, MeshP, Control, Solution)
% Plot solutions at each time step u-p formulation

% initialize figure
figure(1);
tiledlayout(1,3);

% solid displacement
nexttile
plot(MeshU.coords(Control.ploturow), Solution.u(Control.ploturow), 'm', 'LineWidth', 1.5);
hold on
grid on
title('Solid displacement');
xlabel ('coord. [m]');
ylabel('us [m]');
hold off

% solid velocity
nexttile
plot(MeshU.coords(Control.ploturow), Solution.udot(Control.ploturow), 'b', 'LineWidth', 1.5);
hold on
grid on
title('Solid velocity');
xlabel ('coord. [m]');
ylabel('us_dot [m/s]');
hold off

% fluid pressure
nexttile
plot(MeshP.coords(Control.plotprow), Solution.p(Control.plotprow), 'g', 'LineWidth', 1.5);
hold on
grid on
title('Pressure');
xlabel ('coord. [m]');
ylabel('pressure');
hold off

end