% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function [SolutionFreq] = SolverFreqDyn_UP(phi_u, omega2_u, phi_p, omega2_p, Kuu, Kup, Kpp, M, Mhat, S, fu, fp, BC, Control, Iteration)
% Solve linear system for dynamic case

%% Iteration data
uF_old = Iteration.uF_old;
pF_old = Iteration.pF_old;
udot_old = Iteration.uFdot_old;
u2dot_old = Iteration.uF2dot_old;
pdot_old = Iteration.pFdot_old;

xuF_old = Iteration.xuF_old;
xpF_old = Iteration.xpF_old;
xuFdot_old = Iteration.xuFdot_old;
xuF2dot_old = Iteration.xuF2dot_old;
xpFdot_old = Iteration.xpFdot_old;

% time step
dt = Control.dt;

% time integration parameters
beta = Control.beta;
gamma = Control.gamma;
theta = Control.theta;

% global matrix partitioning
KuuFF = Kuu(BC.free_u, BC.free_u);
KupFF = Kup(BC.free_u, BC.free_p);
KpuFF = Kup(BC.free_u, BC.free_p).';
KppFF = Kpp(BC.free_p, BC.free_p);
MFF = M(BC.free_u, BC.free_u);
SFF = S(BC.free_p, BC.free_p);
MhatFF = Mhat(BC.free_p, BC.free_u);

%% First step
if Control.step == 1
    u2dot_old(BC.free_u) = MFF\(-KuuFF*uF_old(BC.free_u) + KupFF*pF_old(BC.free_p) + fu(BC.free_u));
    pdot_old(BC.free_p) = SFF\(MhatFF*u2dot_old(BC.free_u) - KpuFF*udot_old(BC.free_u) -KppFF*pF_old(BC.free_p) + fp(BC.free_p));
    
    xuF2dot_old = ((phi_u).' * MFF * phi_u) \ ((phi_u).' * fu(BC.free_u) - (phi_u).' * KuuFF * phi_u * xuF_old +...
        (phi_u).' * KupFF * phi_p * xpF_old);
    xpFdot_old = ((phi_p).' * SFF * phi_p) \ ((phi_p).' * fp(BC.free_p) + (phi_p).' * MhatFF * phi_u * xuF2dot_old - ...
        (phi_p).' * KpuFF * phi_u * xuFdot_old - (phi_p).' * KppFF * phi_p * xpF_old);
end
%% Matrix partitioning
% matrices time discretization
% Kuubar = (phi_u.') * MFF * phi_u * 4/dt^2 + (phi_u.') * KuuFF * phi_u;
Kuubar = (phi_u.') * MFF * phi_u ./(beta*dt^2) + (phi_u.') * KuuFF * phi_u;
Kupbar = -(phi_u.') * KupFF * phi_p;
% Kpubar = -(phi_p.') * MhatFF * phi_u * 4/dt^2 + (phi_p.') * KpuFF * phi_u * 2/dt;
Kpubar = -(phi_p.') * MhatFF * phi_u ./(beta*dt^2) + (phi_p.') * KpuFF * phi_u * gamma/(beta*dt);
% Kppbar = (phi_p.') * KppFF * phi_p + (phi_p.') * SFF * phi_p * 2/dt;
Kppbar = (phi_p.') * KppFF * phi_p + (phi_p.') * SFF * phi_p ./(theta*dt);

test1 = (phi_u.') * MFF * phi_u;
test2 = (phi_u.') * KuuFF * phi_u;
test3 = (phi_p.') * KppFF * phi_p;

% auxiliar terms for external forces vector
% fubar = (phi_u.') * fu(BC.free_u) + (phi_u.') * MFF * phi_u * (xuF_old * 4/dt^2 + xuFdot_old * 4/dt + xuF2dot_old);
% fpbar = (phi_p.') * fp(BC.free_p) - (phi_p.') * MhatFF * phi_u * (xuF_old * 4/dt^2 + xuFdot_old * 4/dt + xuF2dot_old) + ...
%     (phi_p.') * KpuFF * phi_u * (xuFdot_old + xuF_old * 2/dt) + (phi_p.') * SFF * phi_p * (xpFdot_old + xpF_old * 2/dt);

fubar = (phi_u.') * fu(BC.free_u) + (phi_u.') * MFF * phi_u * (xuF_old ./(beta*dt^2) + xuFdot_old ./(beta*dt) + (1/(2*beta) -1) * xuF2dot_old);
fpbar = (phi_p.') * fp(BC.free_p) - (phi_p.') * MhatFF * phi_u * (xuF_old ./(beta*dt^2) + xuFdot_old ./(beta*dt) + (1/(2*beta)-1) * xuF2dot_old) + ...
    (phi_p.') * KpuFF * phi_u * (xuFdot_old * (gamma/beta -1) + xuF_old * gamma/(beta*dt)) + (phi_p.') * SFF * phi_p * ((1/theta -1) * xpFdot_old + xpF_old ./(theta *dt));

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
uE = BC.fixed_u_value(Control.t);
pE = BC.fixed_p_value(Control.t);

% store u/p fields
u(BC.fixed_u, 1) = uE;
u(BC.free_u, 1) = uF;
p(BC.fixed_p, 1) = pE;
p(BC.free_p, 1) = pF;

%% Gradients
% xuFdot = - xuFdot_old + (xuF - xuF_old) * 2/dt;
% xpFdot = - xpFdot_old + (xpF - xpF_old) * 2/dt;
% xuF2dot = xuF * 4/dt^2 - xuF_old * 4/dt^2 - xuFdot_old * 4/dt - xuF2dot_old;
% 
% udot = - udot_old + (u - uF_old) * 2/dt;
% pdot = - pdot_old + (p - pF_old) * 2/dt;
% u2dot = u * 4/dt^2 - uF_old * 4/dt^2 - udot_old * 4/dt - u2dot_old;

xuFdot = - xuFdot_old * (gamma/beta -1) + (xuF - xuF_old) *gamma/(beta*dt) - xuF2dot_old * dt * (gamma/(2*beta)-1);
xpFdot = - xpFdot_old * (1/theta - 1) + (xpF - xpF_old) ./(theta*dt);
xuF2dot = (xuF - xuF_old) ./(beta*dt^2) - xuFdot_old ./(beta*dt) - xuF2dot_old * (1/(2*beta) -1);

udot = - udot_old * (gamma/beta -1)  + (u - uF_old) * gamma/(beta*dt) - u2dot_old * dt * (gamma/(2*beta)-1);
pdot = - pdot_old * (1/theta - 1) + (p - pF_old) ./(theta*dt);
u2dot = (u - uF_old) ./(beta*dt^2) - udot_old ./(beta*dt) - u2dot_old * (1/(2*beta) -1);

%% Store variables
SolutionFreq.xuF = xuF;
SolutionFreq.xpF = xpF;
SolutionFreq.xuFdot = xuFdot;
SolutionFreq.xpFdot = xpFdot;
SolutionFreq.xuF2dot = xuF2dot;

SolutionFreq.uF = u;
SolutionFreq.uFdot = udot;
SolutionFreq.pF = p;
SolutionFreq.pFdot = pdot;
SolutionFreq.uF2dot = u2dot;

end