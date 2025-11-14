function getAmplificationMatrix(Material, MeshU, MeshP, QuadU, QuadP, Control)
% Get amplification matrix for spectral analysis of time integration scheme

fprintf('Get amplification matrix...');

% time step range
dtmin = 1e-5;
dtmax = 1e-4;
dt = linspace(dtmin, dtmax, 10);
% source frequency
f = 1e3;
% initialize sigma vector
sigma_max = zeros(size(dt));

% time integration parameters
beta = Control.beta;
gamma = Control.gamma;
theta = Control.theta;
alpha = theta;

% [Mss, Msf, Kss, Ksp, Mfs, Mff, Kfp, Kff, Cps, Kpf, Cpp, Cfs] = ComputeMatricesDyn_dCS_UPW_dimensionless(Material, MeshU, MeshP, QuadU, QuadP);
[Mss, Msf, Kss, Ksp, Mfs, Mff, Kfp, Kff, Cps, Kpf, Cpp, Cfs] = ComputeMatricesDyn_BT_UPW(Material, MeshU, MeshP, QuadU, QuadP);

% loop over time steps
for i = 1:length(dt)
    % compute matrices in G = A^-1 * B
    [A, B] = getABmatrices(Mss, Msf, Kss, Ksp, Mfs, Mff, Kfp, Kff, Cps, Kpf, Cpp, Cfs, beta, gamma, theta, alpha, dt(i));

    % amplification matrix
    G = A\B;
    % G = MatrixInvert(A, B, 6);
    % symmetric matrix
    % Gsym = G*G';
    k = 10; % a few largest singular values
    sigma = svds(G, k);
    sigma_max(i) = max(sigma);

    % n = size(B,2);   % number of columns of B = size of input state
    % k = 5;           % number of largest singular values you want
    % Gfun = @(x) A \ (B*x);  % forward operator
    % sigma = svds(Gfun, n, k, 'largest');
    % sigma_max = max(sigma);

    % get eigenvalues
    % eigG = eig(Gsym);
    % eig_max = eigs(Gsym, 1);
    % sigma_max(i) = sqrt(real(eig_max));
end

figure;
plot(dt*f, sigma_max, 'LineWidth', 1.5);
xlabel('\Delta t / T');
ylabel('max singular value \sigma_{max}(G)');
yline(1, '--r', '|\sigma|=1');
title('Maximum one-step amplification');
grid on

end


function [A, B] = getABmatrices(Mss, Msf, Kss, Ksp, Mfs, Mff, Kfp, Kff, Cps, Kpf, Cpp, Cfs, beta, gamma, theta, alpha, dt)
% get A, B matrix from the system
% A*d^(n+1) = F^(n+1) + B*d^(n)
% where d = [us, us_dot, w, p, us_2dot, w_dot, p_dot]

% sizes
ns = length(Mss);
nf = length(Mff);
np = length(Cpp);

% identity matrices
Is = speye(ns);
If = speye(nf);
Ip = speye(np);

% zero matrices
Zss = sparse(ns, ns);
Zff = sparse(nf, nf);
Zpp = sparse(np, np);
Zsf = Zss; Zfs = Zss;
Zsp = sparse(ns, np);
Zps = sparse(np, ns);
Zfp = sparse(nf, np);
Zpf = sparse(np, nf);

% matrix A (7x7)
A = [Is, Zss, Zsf, Zsp, -dt*beta*Is, Zsf, Zsp;
    Zss, Is, Zsf, Zsp, -dt*gamma*Is, Zsf, Zsp;
    Zfs, Zfs, If, Zfp, Zfs, -dt*theta*If, Zfp;
    Zps, Zps, Zpf, Ip, Zps, Zpf, -dt*alpha*Ip;
    Kss, Zss, Zsf, -Ksp, Mss, Msf, Zsp;
    Zfs, Cfs, Kff, -Kfp, Mfs, Mff, Zfp;
    Zps, Cps, Kpf, Zpp, Zps, Zpf, Cpp];

% matrix B (7x7)
B = [Is, dt*Is, Zsf, Zsp, dt^2*(0.5-beta)*Is, Zsf, Zsp;
    Zss, Is, Zsf, Zsp, dt*(1-gamma)*Is, Zsf, Zsp;
    Zfs, Zfs, If, Zfp, Zfs, dt*(1-theta)*If, Zfp;
    Zps, Zps, Zpf, Ip, Zps, Zpf, dt*(1-alpha)*Ip;
    Zss, Zss, Zsf, Zsp, Zss, Zsf, Zsp;
    Zfs, Zfs, Zff, Zfp, Zfs, Zff, Zfp;
    Zps, Zps, Zpf, Zpp, Zps, Zpf, Zpp];
end