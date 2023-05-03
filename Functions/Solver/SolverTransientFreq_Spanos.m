function [SolutionFreq] = SolverTransientFreq_Spanos(phi_u, omega2_u, phi_p, omega2_p, phi_n, omega2_n, Kuu, Kup, Kpu, Kpp, S, Kpn, Knu, Knp, Knn, fu, fp, BC, Control, Iteration)
% ------------------------------------------------------------------------
% Solve linear system for quasi-steady case
% ------------------------------------------------------------------------
% Input parameters: coupled matrices, BC, Control, Iteration
% ------------------------------------------------------------------------
% Outputs:
%   u: solid displacement
%   udot: solid velocity
%   p: fluid pressure
% ------------------------------------------------------------------------

%% Iteration data
uF_old = Iteration.uF_old;
pF_old = Iteration.pF_old;
nF_old = Iteration.nF_old;
xuF_old = Iteration.xuF_old;
xpF_old = Iteration.xpF_old;
xnF_old = Iteration.xnF_old;

% time step
dt = Control.dt;

%% Matrix partitioning
% matrices time discretization
Kuubar = (phi_u.') * Kuu(BC.free_u, BC.free_u) * phi_u;
Kupbar = -(phi_u.') * Kup(BC.free_u, BC.free_p) * phi_p;
Kunbar = sparse(length(BC.free_u), length(BC.free_n));

Kpubar = (phi_p.') * Kpu(BC.free_p, BC.free_u) * phi_u;
Kppbar = (phi_p.') * Kpp(BC.free_p, BC.free_p) * phi_p * dt + (phi_p.') * S(BC.free_p, BC.free_p) * phi_p;
Kpnbar = (phi_p.') * Kpn(BC.free_p, BC.free_n) * phi_n;

Knubar = (phi_n.') * Knu(BC.free_n, BC.free_u) * phi_u;
Knpbar = (phi_n.') * Knp(BC.free_n, BC.free_p) * phi_p * dt;
Knnbar = (phi_n.') * Knn(BC.free_n, BC.free_n);

% auxiliar terms for external forces vector
fubar = (phi_u.') * fu(BC.free_u);
fpbar = (phi_p.') * fp(BC.free_p) * dt + (phi_p.') * Kpu(BC.free_p, BC.free_u) * phi_u * xuF_old +...
    (phi_p.') * S(BC.free_p, BC.free_p) * phi_p * xpF_old + (phi_p.') * Kpn(BC.free_p, BC.free_n) * phi_n * xnF_old;
fnbar = (phi_n.') * Knu(BC.free_n, BC.free_u) * phi_u * xuF_old + (phi_n.') * Knn(BC.free_n, BC.free_n) * phi_n * xnF_old;

% matrices of the system to solve
AFF = [Kuubar, Kupbar, Kunbar;
    Kpubar, Kppbar, Kpnbar;
    Knubar, Knpbar, Knnbar];
BF = [fubar; fpbar; fnbar];
     
%% Solve linear system
% solve for displacement and pressure
x = AFF\BF;

%% Separate u/p
% unknown displacement
xuF = x(1:length(BC.free_u), 1);
% unknown pressure
xpF = x(length(BC.free_u)+1 : length(BC.free_u) + length(BC.free_p), 1);
% unknown porosity
xnF = x(length(BC.free_u) + length(BC.free_p)+1 : end, 1);

%% Mode superposition
% return to original unkowns
uF = phi_u * xuF;
pF = phi_p * xpF;
nF = phi_n * xnF;

% fixed DOFs
uE = BC.fixed_u_value;
pE = BC.fixed_p_value;
nE = BC.fixed_n_value;

% store u/p fields
u(BC.fixed_u, 1) = uE;
u(BC.free_u, 1) = uF;
p(BC.fixed_p, 1) = pE;
p(BC.free_p, 1) = pF;
n(BC.fixed_n, 1) = nE;
n(BC.free_n, 1) = nF;

%% Gradients
xuFdot = (xuF - xuF_old)./dt;
xpFdot = (xpF - xpF_old)./dt;
xnFdot = (xnF - xnF_old)./dt;

udot = (u - uF_old)./dt;
pdot = (p - pF_old)./dt;
ndot = (n - nF_old)./dt;

%% Store variables
SolutionFreq.xuF = xuF;
SolutionFreq.xpF = xpF;
SolutionFreq.xuFdot = xuFdot;
SolutionFreq.xpFdot = xpFdot;
SolutionFreq.xnF = xnF;
SolutionFreq.xnFdot = xnFdot;

SolutionFreq.uF = u;
SolutionFreq.uFdot = udot;
SolutionFreq.pF = p;
SolutionFreq.pFdot = pdot;
SolutionFreq.nF = n;
SolutionFreq.nFdot = ndot;

end