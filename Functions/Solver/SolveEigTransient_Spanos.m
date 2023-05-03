function [phi_ucoupl, omega2_ucoupl, phi_pcoupl, omega2_pcoupl, phi_ncoupl, omega2_ncoupl] = SolveEigTransient_Spanos(Kuu, Kup, Kpp, Knp, Kpu, S, Kpn, Knu, Knn, MeshU, MeshP, MeshN, BC, Control)
% solve eigenproblem for transient systems
% ------------------------------------------------------------------------

% partitioned matrices for unknown variables
KuuFF = Kuu(BC.free_u, BC.free_u);
KppFF = Kpp(BC.free_p, BC.free_p);
KupFF = Kup(BC.free_u, BC.free_p);
KnpFF = Knp(BC.free_n, BC.free_p);

KpuFF = Kpu(BC.free_p, BC.free_u);
SFF = S(BC.free_p, BC.free_p);
KpnFF = Kpn(BC.free_p, BC.free_n);
KnuFF = Knu(BC.free_n, BC.free_u);
KnnFF = Knn(BC.free_n, BC.free_n);

%% Coupled problem
% coupled matrix
K = [KuuFF, -KupFF, sparse(length(BC.free_u), length(BC.free_n));
    sparse(length(BC.free_p), length(BC.free_u)), KppFF, sparse(length(BC.free_p), length(BC.free_n));
    sparse(length(BC.free_n), length(BC.free_u)), KnpFF, sparse(length(BC.free_n), length(BC.free_n))];
C = [sparse(length(BC.free_u), length(BC.free_u)), sparse(length(BC.free_u), length(BC.free_p)), sparse(length(BC.free_u), length(BC.free_n));
    KpuFF, SFF, KpnFF;
    KnuFF, sparse(length(BC.free_n), length(BC.free_p)), KnnFF];
% coupled solution
[phi_coupled, omega2_coupled] = eigenshuffle(full(K));
omega2_coupled = diag(omega2_coupled);

% coupled damped solution
% [phidamp_coupled, omegadamp_coupled] = eig(full(K), full(C));
[phidamp_coupled, omegadamp_coupled] = eigenshuffle(full(K), full(C));
omegadamp_coupled = diag(omegadamp_coupled);

%% Plots coupled
% coupled modes selection
% displacement
phi_ucoupl = phi_coupled(1:length(BC.free_u), 1:length(BC.free_u));
omega2_ucoupl = omega2_coupled(1:length(BC.free_u), 1:length(BC.free_u));
% pressure
phi_pcoupl = phi_coupled(length(BC.free_u)+1:length(BC.free_u)+length(BC.free_p), length(BC.free_u)+1:length(BC.free_u)+length(BC.free_p));
omega2_pcoupl = omega2_coupled(length(BC.free_u)+1:length(BC.free_u)+length(BC.free_p), length(BC.free_u)+1:length(BC.free_u)+length(BC.free_p));
% porosity
phi_ncoupl = phi_coupled(length(BC.free_u)+length(BC.free_p)+1:end, length(BC.free_u)+length(BC.free_p)+1:end);
omega2_ncoupl = omega2_coupled(length(BC.free_u)+length(BC.free_p)+1:end, length(BC.free_u)+length(BC.free_p)+1:end);


% % displacement
% phi_ncoupl = phi_coupled(1:length(BC.free_n), 1:length(BC.free_n));
% omega2_ncoupl = omega2_coupled(1:length(BC.free_n), 1:length(BC.free_n));
% % pressure
% phi_ucoupl = phi_coupled(length(BC.free_n)+1:length(BC.free_n)+length(BC.free_u), length(BC.free_n)+1:length(BC.free_n)+length(BC.free_u));
% omega2_ucoupl = omega2_coupled(length(BC.free_n)+1:length(BC.free_n)+length(BC.free_u), length(BC.free_n)+1:length(BC.free_n)+length(BC.free_u));
% % porosity
% phi_pcoupl = phi_coupled(length(BC.free_n)+length(BC.free_u)+1:end, length(BC.free_n)+length(BC.free_u)+1:end);
% omega2_pcoupl = omega2_coupled(length(BC.free_n)+length(BC.free_u)+1:end, length(BC.free_n)+length(BC.free_u)+1:end);


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
% sort modes porosity
omega2_ncoupl_vec = diag(omega2_ncoupl);
[omega2_ncoupl_sorted, indexcoupl] = sort(omega2_ncoupl_vec);
phi_ncoupl_sorted = zeros(length(phi_ncoupl), length(phi_ncoupl));
for i = 1: length(phi_ncoupl)
    phi_ncoupl_sorted(:,i) = phi_ncoupl(:,indexcoupl(i));
end

% plot together
figure;
hold on
sgtitle(sprintf('PM Mode Shapes at t = %.2f s - COUPLED EIGENPROBLEM', Control.tend));
modeshapeu = zeros(MeshU.nDOF,6);
modeshapep = zeros(MeshP.nDOF,6);
modeshapen = zeros(MeshN.nDOF,6);
for mode = 1:6
    modeshapeu(BC.free_u, mode) = phi_ucoupl_sorted(:,mode);
    modeshapep(BC.free_p, mode) = phi_pcoupl_sorted(:,mode);
    modeshapen(BC.free_n, mode) = phi_ncoupl_sorted(:,mode);
    subplot(2,3, mode);
    plot(MeshU.coords, modeshapeu(:,mode), 'b', 'LineWidth', 2);
    hold on
    plot(MeshP.coords, modeshapep(:,mode), 'k', 'LineWidth', 2);
    plot(MeshN.coords, modeshapen(:,mode), 'r', 'LineWidth', 2);
    xlabel('Length (m)');
    ylabel('Shape');
    legend('u', 'p', 'n');
    title(sprintf('# %.0f, omU = %.2f Hz, omP = %.2d Hz, omN = %.2d Hz', mode, sqrt(omega2_ucoupl_sorted(mode)), sqrt(omega2_pcoupl_sorted(mode)*1e-9), sqrt(omega2_ncoupl_sorted(mode))));
end
hold off

%% Plots dampled coupled
% coupled damped modes selection
% displacement
phidamp_ucoupl = phidamp_coupled(1:length(BC.free_u), 1:length(BC.free_u));
omegadamp_ucoupl = omegadamp_coupled(1:length(BC.free_u), 1:length(BC.free_u));
omega2damp_ucoupl = omegadamp_ucoupl.^2;
% pressure
phidamp_pcoupl = phidamp_coupled(length(BC.free_u)+1:length(BC.free_u)+length(BC.free_p), length(BC.free_u)+1:length(BC.free_u)+length(BC.free_p));
omegadamp_pcoupl = omegadamp_coupled(length(BC.free_u)+1:length(BC.free_u)+length(BC.free_p), length(BC.free_u)+1:length(BC.free_u)+length(BC.free_p));
omega2damp_pcoupl = omegadamp_pcoupl.^2;
% porosity
phidamp_ncoupl = phidamp_coupled(length(BC.free_u)+length(BC.free_p)+1:end, length(BC.free_u)+length(BC.free_p)+1:end);
omegadamp_ncoupl = omegadamp_coupled(length(BC.free_u)+length(BC.free_p)+1:end, length(BC.free_u)+length(BC.free_p)+1:end);
omega2damp_ncoupl = omegadamp_ncoupl.^2;

% sort modes displacement
omega2damp_ucoupl_vec = diag(omega2damp_ucoupl);
[omega2damp_ucoupl_sorted, indexdampcoupl] = sort(omega2damp_ucoupl_vec);
phidamp_ucoupl_sorted = zeros(length(phidamp_ucoupl), length(phidamp_ucoupl));
for i = 1: length(phidamp_ucoupl)
    phidamp_ucoupl_sorted(:,i) = phidamp_ucoupl(:,indexdampcoupl(i));
end
% sort modes pressure
omega2damp_pcoupl_vec = diag(omega2damp_pcoupl);
[omega2damp_pcoupl_sorted, indexdampcoupl] = sort(omega2damp_pcoupl_vec);
phidamp_pcoupl_sorted = zeros(length(phidamp_pcoupl), length(phidamp_pcoupl));
for i = 1: length(phidamp_pcoupl)
    phidamp_pcoupl_sorted(:,i) = phidamp_pcoupl(:,indexdampcoupl(i));
end
% sort modes porosity
omega2damp_ncoupl_vec = diag(omega2damp_ncoupl);
[omega2damp_ncoupl_sorted, indexdampcoupl] = sort(omega2damp_ncoupl_vec);
phidamp_ncoupl_sorted = zeros(length(phidamp_ncoupl), length(phidamp_ncoupl));
for i = 1: length(phidamp_ncoupl)
    phidamp_ncoupl_sorted(:,i) = phidamp_ncoupl(:,indexdampcoupl(i));
end

% plot together
figure;
hold on
sgtitle(sprintf('PM Mode Shapes at t = %.2f s - DAMPED COUPLED EIGENPROBLEM', Control.tend));
modeshapeu = zeros(MeshU.nDOF,6);
modeshapep = zeros(MeshP.nDOF,6);
modeshapen = zeros(MeshN.nDOF,6);
for mode = 1:6
    modeshapeu(BC.free_u, mode) = phidamp_ucoupl_sorted(:,mode);
    modeshapep(BC.free_p, mode) = phidamp_pcoupl_sorted(:,mode);
    modeshapen(BC.free_n, mode) = phidamp_ncoupl_sorted(:,mode);
    subplot(2,3, mode);
    plot(MeshU.coords, modeshapeu(:,mode), 'b', 'LineWidth', 2);
    hold on
    plot(MeshP.coords, modeshapep(:,mode), 'k', 'LineWidth', 2);
    plot(MeshN.coords, modeshapen(:,mode), 'r', 'LineWidth', 2);
    xlabel('Length (m)');
    ylabel('Shape');
    legend('u', 'p', 'n');
    title(sprintf('# %.0f, omU = %.2f Hz, omP = %.2d Hz, omN = %.2d Hz', mode, sqrt(omega2damp_ucoupl_sorted(mode)), sqrt(omega2damp_pcoupl_sorted(mode)*1e-9), sqrt(omega2damp_ncoupl_sorted(mode))));
end
hold off

end