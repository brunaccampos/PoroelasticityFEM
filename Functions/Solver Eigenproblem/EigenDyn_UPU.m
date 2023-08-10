function [phi_u, omega2_u, phi_p, omega2_p, phi_uf, omega2_uf] = EigenDyn_UPU(Kss, Ksp, Mss, Kpf, Kps, Kpp, Kfp, Mff, Msf, Mfs, MeshU, MeshP, BC, Control)
% solve eigenproblem for dynamic systems considering Biot theory
% ------------------------------------------------------------------------

% partitioned matrices for unknown variables
MssFF = Mss(BC.free_u, BC.free_u);
MsfFF = Msf(BC.free_u, BC.free_u);
MfsFF = Mfs(BC.free_u, BC.free_u);
MffFF = Mff(BC.free_u, BC.free_u);

KssFF = Kss(BC.free_u, BC.free_u);
KspFF = Ksp(BC.free_u, BC.free_p);
KpsFF = Kps(BC.free_p, BC.free_u);
KppFF = Kpp(BC.free_p, BC.free_p);
KpfFF = Kpf(BC.free_p, BC.free_u);
KfpFF = Kfp(BC.free_u, BC.free_p);

ZuuFF = sparse(length(BC.free_u), length(BC.free_u));
ZupFF = sparse(length(BC.free_u), length(BC.free_p));
ZpuFF = ZupFF.';
ZppFF = sparse(length(BC.free_p), length(BC.free_p));

%% Uncoupled problem
% solid displacement uncoupled
[phi_u, omega2_u] = eig(full(KssFF), full(MssFF));
% pressure uncoupled
[phi_p, omega2_p] = eig(full(KppFF));

%% Coupled problem
% coupled matrix
Mbar = [MssFF, ZupFF, MsfFF;
    ZpuFF, ZppFF, ZpuFF;
    MfsFF, ZupFF, MffFF];
Kbar = [KssFF, -KspFF, ZuuFF;
    KpsFF, KppFF, KpfFF;
    ZuuFF, -KfpFF, ZuuFF];
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

%% Plots coupled
% coupled modes selection
% solid displacement
phi_ucoupl = phi_coupled(1:length(BC.free_u), 1:length(BC.free_u));
omega2_ucoupl = omega2_coupled(1:length(BC.free_u), 1:length(BC.free_u));
% pressure
phi_pcoupl = phi_coupled(length(BC.free_u)+1:length(BC.free_u)+length(BC.free_p), length(BC.free_u)+1:length(BC.free_u)+length(BC.free_p));
omega2_pcoupl = omega2_coupled(length(BC.free_u)+1:length(BC.free_u)+length(BC.free_p), length(BC.free_u)+1:length(BC.free_u)+length(BC.free_p));
% fluid displacement
phi_ufcoupl = phi_coupled(length(BC.free_u)+length(BC.free_p)+1:end, length(BC.free_u)+length(BC.free_p)+1:end);
omega2_ufcoupl = omega2_coupled(length(BC.free_u)+length(BC.free_p)+1:end, length(BC.free_u)+length(BC.free_p)+1:end);

% sort modes solid displacement
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
% sort modes fluid displacement
omega2_ufcoupl_vec = diag(omega2_ufcoupl);
[omega2_ufcoupl_sorted, indexcoupl] = sort(omega2_ufcoupl_vec);
phi_ufcoupl_sorted = zeros(length(phi_ufcoupl), length(phi_ufcoupl));
for i = 1: length(phi_ufcoupl)
    phi_ufcoupl_sorted(:,i) = phi_ufcoupl(:,indexcoupl(i));
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
    legend('u_s', 'p');
    title(sprintf('# %.0f, omU = %.2f Hz, omP = %.2d Hz', mode, sqrt(omega2_ucoupl_sorted(mode)), sqrt(omega2_pcoupl_sorted(mode)*1e-9)));
end
hold off

end