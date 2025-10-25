% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function [phi_ucoupl, omega2_ucoupl, phi_pcoupl, omega2_pcoupl] = EigenTr_UP(Kuu, Kup, Kpp, S, MeshU, MeshP, BC, Control)
% solve eigenproblem for transient systems

% partitioned matrices for unknown variables
KuuFF = Kuu(BC.free_u, BC.free_u);
KppFF = Kpp(BC.free_p, BC.free_p);
KupFF = Kup(BC.free_u, BC.free_p);
KpuFF = KupFF.';
SFF = S(BC.free_p, BC.free_p);

%% Uncoupled problem
% displacement uncoupled
[phi_u, omega2_u] = eig(full(KuuFF));
% pressure uncoupled
[phi_p, omega2_p] = eig(full(KppFF));

%% Coupled problem
% coupled matrix
K = [KuuFF, -KupFF;
    sparse(length(KppFF), length(KuuFF)), KppFF];
C = [sparse(length(KuuFF), length(KuuFF)), sparse(length(KuuFF), length(KppFF));
    KpuFF, SFF];
% coupled solution
[phi_coupled, omega2_coupled] = eig(full(K));
% coupled damped solution
[phidamp_coupled, omegadamp_coupled] = eig(full(K), full(C));

%% Plots uncoupled
% sort modes displacement
omega2_u_vec = diag(omega2_u);
[omega2_u_sorted, index] = sort(omega2_u_vec);
phi_u_sorted = zeros(length(phi_u), length(phi_u));
for i = 1: length(phi_u)
    phi_u_sorted(:,i) = phi_u(:,index(i));
end
% sort modes pressure
omega2_p_vec = diag(omega2_p);
[omega2_p_sorted, index] = sort(omega2_p_vec);
phi_p_sorted = zeros(length(phi_p), length(phi_p));
for i = 1: length(phi_p)
    phi_p_sorted(:,i) = phi_p(:,index(i));
end
% plot together
figure;
hold on
sgtitle(sprintf('PM Mode Shapes at t = %.2f s - UNCOUPLED EIGENPROBLEM', Control.tend));
modeshapeu = zeros(MeshU.nDOF,6);
modeshapep = zeros(MeshP.nDOF,6);
for mode = 1:6
    modeshapeu(BC.free_u, mode) = phi_u_sorted(:,mode);
    modeshapep(BC.free_p, mode) = phi_p_sorted(:,mode);
    subplot(2,3, mode);
    plot(MeshU.coords, modeshapeu(:,mode), 'b', 'LineWidth', 2);
    hold on
    plot(MeshP.coords, modeshapep(:,mode), 'k', 'LineWidth', 2);
    xlabel('Length (m)');
    ylabel('Shape');
    legend('u', 'p');
    title(sprintf('# %.0f, omU = %.2f Hz, omP = %.2d Hz', mode, sqrt(omega2_u_sorted(mode)*1e9), sqrt(omega2_p_sorted(mode)*1e-9)));
end
hold off

%% Plots coupled
% coupled modes selection
phi_ucoupl = phi_coupled(1:length(BC.free_u), 1:length(BC.free_u));
omega2_ucoupl = omega2_coupled(1:length(BC.free_u), 1:length(BC.free_u));
phi_pcoupl = phi_coupled(length(BC.free_u)+1:end, length(BC.free_u)+1:end);
omega2_pcoupl = omega2_coupled(length(BC.free_u)+1:end, length(BC.free_u)+1:end);

% sort modes displacement
omega2_ucoupl_vec = diag(omega2_ucoupl);
[omega2_ucoupl_sorted, indexcoupl] = sort(omega2_ucoupl_vec);
phi_ucoupl_sorted = zeros(length(phi_ucoupl), length(phi_ucoupl));
for i = 1: length(phi_ucoupl)
    phi_ucoupl_sorted(:,i) = phi_ucoupl(:,indexcoupl(i));
end
% sort modes pressure
omega2_pcoupl_vec = diag(omega2_pcoupl);
[omega2_pcoupl_sorted, indexcoupl] = sort(omega2_pcoupl_vec);
phi_pcoupl_sorted = zeros(length(phi_pcoupl), length(phi_pcoupl));
for i = 1: length(phi_pcoupl)
    phi_pcoupl_sorted(:,i) = phi_pcoupl(:,indexcoupl(i));
end
% plot together
figure;
hold on
sgtitle(sprintf('PM Mode Shapes at t = %.2f s - COUPLED EIGENPROBLEM', Control.tend));
modeshapeu = zeros(MeshU.nDOF,6);
modeshapep = zeros(MeshP.nDOF,6);
for mode = 1:6
    modeshapeu(BC.free_u, mode) = phi_ucoupl_sorted(:,mode);
    modeshapep(BC.free_p, mode) = phi_pcoupl_sorted(:,mode);
    subplot(2,3, mode);
    plot(MeshU.coords, modeshapeu(:,mode), 'b', 'LineWidth', 2);
    hold on
    plot(MeshP.coords, modeshapep(:,mode), 'k', 'LineWidth', 2);
    xlabel('Length (m)');
    ylabel('Shape');
    legend('u', 'p');
    title(sprintf('# %.0f, omU = %.2f Hz, omP = %.2d Hz', mode, sqrt(omega2_ucoupl_sorted(mode)*1e9), sqrt(omega2_pcoupl_sorted(mode)*1e-9)));
end
hold off

%% Plots damped coupled
% damped coupled modes selection
phidamp_ucoupl = phidamp_coupled(1:length(BC.free_u), 1:length(BC.free_u));
omegadamp_ucoupl = omegadamp_coupled(1:length(BC.free_u), 1:length(BC.free_u));
omega2damp_ucoupl = omegadamp_ucoupl.^2;
phidamp_pcoupl = phidamp_coupled(length(BC.free_u)+1:end, length(BC.free_u)+1:end);
omegadamp_pcoupl = omegadamp_coupled(length(BC.free_u)+1:end, length(BC.free_u)+1:end);
omega2damp_pcoupl = omegadamp_pcoupl.^2;

% sort modes displacement
omegadamp_ucoupl_vec = diag(omegadamp_ucoupl);
[omegadamp_ucoupl_sorted, indexdampcoupl] = sort(omegadamp_ucoupl_vec);
phidamp_ucoupl_sorted = zeros(length(phidamp_ucoupl), length(phidamp_ucoupl));
for i = 1: length(phidamp_ucoupl)
    phidamp_ucoupl_sorted(:,i) = phidamp_ucoupl(:,indexdampcoupl(i));
end
% sort modes pressure
omegadamp_pcoupl_vec = diag(omegadamp_pcoupl);
[omegadamp_pcoupl_sorted, indexdampcoupl] = sort(omegadamp_pcoupl_vec);
phidamp_pcoupl_sorted = zeros(length(phidamp_pcoupl), length(phidamp_pcoupl));
for i = 1: length(phidamp_pcoupl)
    phidamp_pcoupl_sorted(:,i) = phidamp_pcoupl(:,indexdampcoupl(i));
end
% plot together
figure;
hold on
sgtitle(sprintf('PM Mode Shapes at t = %.2f s - DAMPED COUPLED EIGENPROBLEM', Control.tend));
modeshapeu = zeros(MeshU.nDOF,6);
modeshapep = zeros(MeshP.nDOF,6);
for mode = 1:6
    modeshapeu(BC.free_u, mode) = phidamp_ucoupl_sorted(:,mode);
    modeshapep(BC.free_p, mode) = phidamp_pcoupl_sorted(:,mode);
    subplot(2,3, mode);
    plot(MeshU.coords, modeshapeu(:,mode), 'b', 'LineWidth', 2);
    hold on
    plot(MeshP.coords, modeshapep(:,mode), 'k', 'LineWidth', 2);
    xlabel('Length (m)');
    ylabel('Shape');
    legend('u', 'p');
    title(sprintf('# %.0f, omU = %.2f Hz, omP = %.2d Hz', mode, omegadamp_ucoupl_sorted(mode)*1e9, omegadamp_pcoupl_sorted(mode)*1e-9));
end
hold off

end