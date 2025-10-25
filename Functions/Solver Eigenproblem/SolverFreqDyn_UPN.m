% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function [SolutionFreq] = SolverFreqDyn_UPN(phi_u, omega2_u, phi_p, omega2_p, phi_n, omega2_n, Kuu, Kup, Kpu, Kpp, S, Kpn, Knu, Knp, Knn, Muu, Mpu, Mnu, fu, fp, BC, Control, Iteration)
% Solve linear system for dynamic case

%% Iteration data
uF_old = Iteration.uF_old;
pF_old = Iteration.pF_old;
nF_old = Iteration.nF_old;
udot_old = Iteration.uFdot_old;
u2dot_old = Iteration.uF2dot_old;
pdot_old = Iteration.pFdot_old;
ndot_old = Iteration.nFdot_old;

xuF_old = Iteration.xuF_old;
xpF_old = Iteration.xpF_old;
xnF_old = Iteration.xnF_old;
xuFdot_old = Iteration.xuFdot_old;
xuF2dot_old = Iteration.xuF2dot_old;
xpFdot_old = Iteration.xpFdot_old;
xnFdot_old = Iteration.xnFdot_old;

% time step
dt = Control.dt;

% time integration parameters
beta = Control.beta;
gamma = Control.gamma;
theta = Control.theta;
lambda = Control.lambda;

% global matrix partitioning
KuuFF = Kuu(BC.free_u, BC.free_u);
KupFF = Kup(BC.free_u, BC.free_p);
MuuFF = Muu(BC.free_u, BC.free_u);

MpuFF = Mpu(BC.free_p, BC.free_u);
KpuFF = Kpu(BC.free_p, BC.free_u);
KppFF = Kpp(BC.free_p, BC.free_p);
SFF = S(BC.free_p, BC.free_p);
KpnFF = Kpn(BC.free_p, BC.free_n);

MnuFF = Mnu(BC.free_n, BC.free_u);
KnuFF = Knu(BC.free_n, BC.free_u);
KnpFF = Knp(BC.free_n, BC.free_p);
KnnFF = Knn(BC.free_n, BC.free_n);

%% First step
if Control.step == 1
    u2dot_old(BC.free_u) = MuuFF\(-KuuFF*uF_old(BC.free_u) + KupFF*pF_old(BC.free_p) + fu(BC.free_u));
    ndot_old(BC.free_n) = KnnFF\(MnuFF*u2dot_old(BC.free_u) - KnuFF*udot_old(BC.free_u) - KnpFF*pF_old(BC.free_p));
    pdot_old(BC.free_p) = SFF\(MpuFF*u2dot_old(BC.free_u) - KpuFF*udot_old(BC.free_u) - KppFF*pF_old(BC.free_p) - KpnFF*ndot_old(BC.free_n) + fp(BC.free_p));
    
    xuF2dot_old = ((phi_u.')*MuuFF*phi_u) \ ((phi_u.')*fu(BC.free_u) - (phi_u.')*KuuFF*phi_u*xuF_old +...
        (phi_u).' * KupFF * phi_p * xpF_old);
    xnFdot_old = ((phi_n.')*KnnFF*phi_n) \ ((phi_n.')*MnuFF*phi_u*xuF2dot_old - (phi_n.')*KnuFF*phi_u*xuFdot_old -...
        (phi_n.')*KnpFF*phi_p*xpF_old);
    xpFdot_old = ((phi_p.')*SFF*phi_p) \ ((phi_p.')*fp(BC.free_p) + (phi_p.')*MpuFF*phi_u*xuF2dot_old - ...
        (phi_p.')*KpuFF*phi_u*xuFdot_old - (phi_p.')*KppFF*phi_p*xpF_old - (phi_p.')*KpnFF*phi_n*xnFdot_old);
end
%% Matrix partitioning
% matrices time discretization
Kuubar = (phi_u.') * MuuFF * phi_u ./(beta*dt^2) + (phi_u.') * KuuFF * phi_u;
Kupbar = -(phi_u.') * KupFF * phi_p;
Kpubar = -(phi_p.') * MpuFF * phi_u ./(beta*dt^2) + (phi_p.') * KpuFF * phi_u * gamma/(beta*dt);
Kppbar = (phi_p.') * KppFF * phi_p + (phi_p.') * SFF * phi_p ./(theta*dt);
Kpnbar = (phi_p.') * KpnFF * phi_n ./(lambda*dt);
Knubar = -(phi_n.') * MnuFF * phi_u + (phi_n.') * KnuFF * phi_u .* (gamma/(beta*dt));
Knpbar = (phi_n.') * KnpFF * phi_p;
Knnbar = (phi_n.') * KnnFF * phi_n ./(lambda*dt);

% auxiliar terms for external forces vector
fubar = (phi_u.') * fu(BC.free_u) + (phi_u.') * MuuFF * phi_u * (xuF_old ./(beta*dt^2) + xuFdot_old ./(beta*dt) + (1/(2*beta) -1) * xuF2dot_old);
fpbar = (phi_p.') * fp(BC.free_p) - (phi_p.') * MpuFF * phi_u * (xuF_old ./(beta*dt^2) + xuFdot_old ./(beta*dt) + (1/(2*beta) -1) * xuF2dot_old) + ...
    (phi_p.') * KpuFF * phi_u * (xuFdot_old * (gamma/beta -1) + xuF_old * gamma/(beta*dt) + dt*(gamma/(2*beta)-1) * xuF2dot_old) +...
    (phi_p.') * SFF * phi_p * ((1/theta -1) * xpFdot_old + xpF_old ./(theta*dt)) + (phi_p.') * KpnFF * phi_n * ((1/lambda-1) *xnFdot_old + xnF_old./(lambda*dt));
fnbar = (phi_n.') * MnuFF * phi_u * (-xuF_old./(beta*dt^2) - xuFdot_old./(beta*dt) + (1-1/(2*beta))*xuF2dot_old) +...
    (phi_n.') * KnuFF * phi_u * ((gamma/beta-1)*xuFdot_old + xuF_old.*gamma/(beta*dt) + xuF2dot_old*dt*(gamma/(2*beta) -1)) +...
    (phi_n.') * KnnFF * phi_n * ((1/lambda -1)*xnFdot_old + xnF_old./(lambda*dt));

% matrices of the system to solve
AFF = [Kuubar, Kupbar, sparse(length(KuuFF), length(KnnFF));
    Kpubar, Kppbar, Kpnbar;
    Knubar, Knpbar, Knnbar];
BF = [fubar; fpbar; fnbar];
     
%% Solve linear system
% solve for displacement and pressure
x = AFF\BF;

%% Separate u/p
% unknown displacement
xuF = x(1:length(BC.free_u),1);
% unknown pressure
xpF = x(length(BC.free_u)+1 : length(BC.free_u) + length(BC.free_p),1);
% unknown porosity
xnF = x(length(BC.free_u)+length(BC.free_p)+1 : end,1);

%% Mode superposition
% return to original unkowns
uF = phi_u * xuF;
pF = phi_p * xpF;
nF = phi_n * xnF;

% fixed DOFs
uE = BC.fixed_u_value(Control.t);
pE = BC.fixed_p_value(Control.t);
nE = BC.fixed_n_value;

% store u/p fields
u(BC.fixed_u, 1) = uE;
u(BC.free_u, 1) = uF;
p(BC.fixed_p, 1) = pE;
p(BC.free_p, 1) = pF;
n(BC.fixed_n, 1) = nE;
n(BC.free_n, 1) = nF;

%% Gradients
xuFdot = - xuFdot_old * (gamma/beta -1) + (xuF - xuF_old) *gamma/(beta*dt) - xuF2dot_old * dt * (gamma/(2*beta)-1);
xpFdot = - xpFdot_old * (1/theta - 1) + (xpF - xpF_old) ./(theta*dt);
xuF2dot = (xuF - xuF_old) ./(beta*dt^2) - xuFdot_old ./(beta*dt) - xuF2dot_old * (1/(2*beta) -1);
xnFdot = (1-1/lambda) * xnFdot_old + (xnF - xnF_old)./(lambda*dt);

udot = - udot_old * (gamma/beta -1)  + (u - uF_old) * gamma/(beta*dt) - u2dot_old * dt * (gamma/(2*beta)-1);
pdot = - pdot_old * (1/theta - 1) + (p - pF_old) ./(theta*dt);
u2dot = (u - uF_old) ./(beta*dt^2) - udot_old ./(beta*dt) - u2dot_old * (1/(2*beta) -1);
ndot = (1-1/lambda) * ndot_old + (n - nF_old)./(lambda*dt);

%% Store variables
SolutionFreq.xuF = xuF;
SolutionFreq.xpF = xpF;
SolutionFreq.xuFdot = xuFdot;
SolutionFreq.xpFdot = xpFdot;
SolutionFreq.xuF2dot = xuF2dot;
SolutionFreq.xnF = xnF;
SolutionFreq.xnFdot = xnFdot;

SolutionFreq.uF = u;
SolutionFreq.uFdot = udot;
SolutionFreq.pF = p;
SolutionFreq.pFdot = pdot;
SolutionFreq.uF2dot = u2dot;
SolutionFreq.nF = n;
SolutionFreq.nFdot = ndot;

end