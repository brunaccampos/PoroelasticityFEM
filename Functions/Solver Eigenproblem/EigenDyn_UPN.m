function [phi_ucoupl, omega2_ucoupl, phi_pcoupl, omega2_pcoupl, phi_ncoupl, omega2_ncoupl] = EigenDyn_UPN(Muu, Mpu, Mnu, Kuu, Kup, Kpp, Knp, MeshU, MeshP, MeshN, BC, Control)
% solve eigenproblem for dynamic systems considering Spanos theory
% ------------------------------------------------------------------------

% partitioned matrices for unknown variables
MuuFF = Muu(BC.free_u, BC.free_u);
MpuFF = Mpu(BC.free_p, BC.free_u);
MnuFF = Mnu(BC.free_n, BC.free_u);

KuuFF = Kuu(BC.free_u, BC.free_u);
KupFF = Kup(BC.free_u, BC.free_p);
KppFF = Kpp(BC.free_p, BC.free_p);
KnpFF = Knp(BC.free_n, BC.free_p);

%% Uncoupled problem
% displacement uncoupled
[phi_u, omega2_u] = eig(full(KuuFF), full(MuuFF));
% pressure uncoupled
[phi_p, omega2_p] = eig(full(KppFF));

%% Coupled problem
% coupled matrix
Mbar = [MuuFF, sparse(length(BC.free_u), length(BC.free_p)), sparse(length(BC.free_u), length(BC.free_n));
    -MpuFF, sparse(length(BC.free_p), length(BC.free_p)), sparse(length(BC.free_p), length(BC.free_n));
    -MnuFF, sparse(length(BC.free_n), length(BC.free_p)), sparse(length(BC.free_n), length(BC.free_n))];

Kbar = [KuuFF, -KupFF, sparse(length(BC.free_u), length(BC.free_n));
    sparse(length(BC.free_p), length(BC.free_u)), KppFF, sparse(length(BC.free_p), length(BC.free_n));
    sparse(length(BC.free_n), length(BC.free_u)), KnpFF, sparse(length(BC.free_n), length(BC.free_n))];

% coupled solution
[phi_coupled, omega2_coupled] = eigenshuffle(full(Kbar), full(Mbar));

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
sgtitle(sprintf('PM Mode Shapes at t = %.2f s - COUPLED PROBLEM', Control.tend));
modeshapeu = zeros(MeshU.nDOF,6);
modeshapep = zeros(MeshP.nDOF,6);
modeshapen = zeros(MeshN.nDOF,6);
for mode = 1:6
    modeshapeu(BC.free_u, mode) = phi_ucoupl_sorted(:,mode) / norm(phi_ucoupl_sorted(:,mode));
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

end