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

% submatrices from A
Ass = Kss + Mss./(beta*dt^2);
Asf = Msf./(theta*dt);
Asp = -Ksp;

Afs = Mfs./(beta*dt^2) + Cfs*gamma/beta/dt;
Aff = Mff./(theta*dt) + Kff;
Afp = -Kfp;

Aps = Cps*gamma/beta/dt;
Apf = Kpf;
App = Cpp./(alpha*dt);

% coefficient matrix in A*d^(n+1)
A = [Ass, Asf, Asp;
    Afs, Aff, Afp;
    Aps, Apf, App];

% submatrices from B
Bs = [Mss*(1/2/beta-1), Mss*(beta/dt), Mss*(beta/dt^2), ...
    zeros(size(Ksp,1), size(Ksp,2)), zeros(size(Ksp,1), size(Ksp,2)), ...
    - Msf*(1-1/theta), Msf*(theta/dt)];

Bf = [Mfs*(1/2/beta-1) + Cfs*dt*(gamma/2/beta-1), Mfs*(beta/dt) + Cfs*(gamma/beta-1), Mfs*(beta/dt^2) + Cfs*(gamma/beta/dt), ...
    zeros(size(Ksp,1), size(Ksp,2)), zeros(size(Ksp,1), size(Ksp,2)), ...
    - Mff*(1-1/theta), Mff*(theta/dt)];

Bp = [Cps*dt*(gamma/2/beta-1), Cps*(gamma/beta-1), Cps*(gamma/beta/dt), ...
    - Cpp*(1-1/alpha), Cpp*(alpha/dt), ...
    zeros(size(Kpf,1), size(Kpf,2)), zeros(size(Kpf,1), size(Kpf,2))];

% coefficient matrix in B*d^(n)
B = [Bs; Bf; Bp];

end