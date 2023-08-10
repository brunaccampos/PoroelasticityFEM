function [phi_u, omega2_u, phi_p, omega2_p] = EigenDyn_UP(Kuu, Kup, Kpp, M, Mhat, MeshU, MeshP, BC, Control)
% solve eigenproblem for dynamic systems considering Biot theory
% ------------------------------------------------------------------------

% partitioned matrices for unknown variables
KuuFF = Kuu(BC.free_u, BC.free_u);
MFF = M(BC.free_u, BC.free_u);
KppFF = Kpp(BC.free_p, BC.free_p);
MhatFF = Mhat(BC.free_p, BC.free_u);
KupFF = Kup(BC.free_u, BC.free_p);

%% Uncoupled problem
% displacement uncoupled
[phi_u, omega2_u] = eig(full(KuuFF), full(MFF));
% pressure uncoupled
[phi_p, omega2_p] = eig(full(KppFF));

%% Coupled problem
% coupled matrix
Mbar = [MFF, sparse(length(BC.free_u), length(BC.free_p));
    -MhatFF, sparse(length(BC.free_p), length(BC.free_p))];
Kbar = [KuuFF, -KupFF;
    sparse(length(BC.free_p), length(BC.free_u)), KppFF];
% coupled solution
[phi_coupled, omega2_coupled] = eig(full(Kbar), full(Mbar));

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
sgtitle(sprintf('PM Mode Shapes at t = %.2f s - UNCOUPLED PROBLEM', Control.tend));
modeshapeu = zeros(MeshU.nDOF,6);
modeshapep = zeros(MeshP.nDOF,6);
for mode = 1:6
    modeshapeu(BC.free_u, mode) = phi_u_sorted(:,mode) / norm(phi_u_sorted(:,mode));
    modeshapep(BC.free_p, mode) = phi_p_sorted(:,mode);
    subplot(2,3, mode);
    plot(MeshU.coords, modeshapeu(:,mode), 'b', 'LineWidth', 2);
    hold on
    plot(MeshP.coords, modeshapep(:,mode), 'k', 'LineWidth', 2);
    xlabel('Length (m)');
    ylabel('Shape');
    legend('u', 'p');
    title(sprintf('# %.0f, omU = %.2f Hz, omP = %.2d Hz', mode, sqrt(omega2_u_sorted(mode)), sqrt(omega2_p_sorted(mode)*1e-9)));
end
hold off

% energy displacement field
figure;
subplot(1,2,1);
plot(omega2_u_sorted, 'b', 'LineWidth', 2);
hold on
xlabel('Mode number');
ylabel('\omega_u^2');
title('Energy of modes: solid displacement');
hold off

% energy pressure field
subplot(1,2,2);
plot(omega2_p_sorted, 'k', 'LineWidth', 2);
hold on
xlabel('Mode number');
ylabel('\omega_p^2');
title('Energy of modes: fluid pressure');
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
sgtitle(sprintf('PM Mode Shapes at t = %.2f s - COUPLED PROBLEM', Control.tend));
modeshapeu = zeros(MeshU.nDOF,6);
modeshapep = zeros(MeshP.nDOF,6);
for mode = 1:6
    modeshapeu(BC.free_u, mode) = phi_ucoupl_sorted(:,mode) / norm(phi_ucoupl_sorted(:,mode));
    modeshapep(BC.free_p, mode) = phi_pcoupl_sorted(:,mode);
    subplot(2,3, mode);
    plot(MeshU.coords, modeshapeu(:,mode), 'b', 'LineWidth', 2);
    hold on
    plot(MeshP.coords, modeshapep(:,mode), 'k', 'LineWidth', 2);
    xlabel('Length (m)');
    ylabel('Shape');
    legend('u', 'p');
    title(sprintf('# %.0f, omU = %.2f Hz, omP = %.2d Hz', mode, sqrt(omega2_ucoupl_sorted(mode)), sqrt(omega2_pcoupl_sorted(mode)*1e-9)));
end
hold off

end