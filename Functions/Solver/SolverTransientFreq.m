function [SolutionFreq] = SolverTransientFreq(phi_u, omega2_u, phi_p, omega2_p, Kuu, Kup, Kpp, S, fu, fp, BC, Control, Iteration)
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
xuF_old = Iteration.xuF_old;
xpF_old = Iteration.xpF_old;

% time step
dt = Control.dt;

%% Matrix partitioning
% matrices time discretization
Kuubar = (phi_u.') * Kuu(BC.free_u, BC.free_u) * phi_u;
Kuubar = omega2_u;
Kupbar = -(phi_u.') * Kup(BC.free_u, BC.free_p) * phi_p;
Kpubar = (phi_p.') * Kup(BC.free_u, BC.free_p)' * phi_u;
Kppbar = (phi_p.') * Kpp(BC.free_p, BC.free_p) * phi_p * dt + (phi_p.') * S(BC.free_p, BC.free_p) * phi_p;
Kppbar = omega2_p * dt + (phi_p.') * S(BC.free_p, BC.free_p) * phi_p;

% auxiliar terms for external forces vector
fubar = (phi_u.') * fu(BC.free_u);
fpbar = (phi_p.') * fp(BC.free_p) * dt + (phi_p.') * Kup(BC.free_u, BC.free_p)' * phi_u * xuF_old +...
    (phi_p.') * S(BC.free_p, BC.free_p) * phi_p * xpF_old;

% matrices of the system to solve
AFF = [Kuubar, Kupbar;
    Kpubar, Kppbar];
BF = [fubar; fpbar];
     
%% Solve linear system
% solve for displacement and pressure
x = AFF\BF;

%% Separate u/p
% unknown displacement
xuF = x(1:length(BC.free_u),1);
% unknown pressure
xpF = x(length(BC.free_u)+1 : length(BC.free_u) + length(BC.free_p),1);

%% Mode superposition
% return to original unkowns
uF = phi_u * xuF;
pF = phi_p * xpF;

% fixed DOFs
uE = BC.fixed_u_value;
pE = BC.fixed_p_value;

% store u/p fields
u(BC.fixed_u, 1) = uE;
u(BC.free_u, 1) = uF;
p(BC.fixed_p, 1) = pE;
p(BC.free_p, 1) = pF;

%% Gradients
xuFdot = (xuF - xuF_old)./dt;
xpFdot = (xpF - xpF_old)./dt;

udot = (u - uF_old)./dt;
pdot = (p - pF_old)./dt;

%% Store variables
SolutionFreq.xuF = xuF;
SolutionFreq.xpF = xpF;
SolutionFreq.xuFdot = xuFdot;
SolutionFreq.xpFdot = xpFdot;

SolutionFreq.uF = u;
SolutionFreq.uFdot = udot;
SolutionFreq.pF = p;
SolutionFreq.pFdot = pdot;

end