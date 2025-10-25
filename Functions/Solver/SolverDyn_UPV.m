% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function [Solution] = SolverDyn_UPV(Mss, Msf, Kss, Ksp, Kps, Kpf, Mff, Kff, Cfs, Kfp, Kpp, fu, fp, ff, BC, Control, Iteration)
% Solve linear system for dynamic case

%% Iteration data
u_old = Iteration.u_old;
udot_old = Iteration.udot_old;
u2dot_old = Iteration.u2dot_old;
p_old = Iteration.p_old;
pdot_old = Iteration.pdot_old;
ufdot_old = Iteration.ufdot_old;
uf2dot_old = Iteration.uf2dot_old;

% current time step
dt = Control.dtc;

% time integration parameters
beta = Control.beta;
gamma = Control.gamma;
theta = Control.theta;
lambda = beta;
alpha = gamma;

%% Matrix partitioning
% matrices time discretization
Kssbar = Kss + Mss./(beta*dt^2);
Kspbar = -Ksp;
Ksfbar = Msf/theta/dt;

Kpsbar = Kps*gamma/beta/dt;
Kpfbar = Kpf;
Kppbar = Kpp*(1/alpha/dt);
Kfsbar = -Cfs*gamma/beta/dt;
Kfpbar = -Kfp;
Kffbar = Mff/theta/dt + Kff;

% matrix partitioning
Kss_EE = Kssbar(BC.fixed_u, BC.fixed_u);
Kss_EF = Kssbar(BC.fixed_u, BC.free_u);
Kss_FE = Kssbar(BC.free_u, BC.fixed_u);
Kss_FF = Kssbar(BC.free_u, BC.free_u);

Ksf_EE = Ksfbar(BC.fixed_u, BC.fixed_ufdot);
Ksf_EF = Ksfbar(BC.fixed_u, BC.free_ufdot);
Ksf_FE = Ksfbar(BC.free_u, BC.fixed_ufdot);
Ksf_FF = Ksfbar(BC.free_u, BC.free_ufdot);

Ksp_EE = Kspbar(BC.fixed_u, BC.fixed_p);
Ksp_EF = Kspbar(BC.fixed_u, BC.free_p);
Ksp_FE = Kspbar(BC.free_u, BC.fixed_p);
Ksp_FF = Kspbar(BC.free_u, BC.free_p);

Kfs_EE = Kfsbar(BC.fixed_ufdot, BC.fixed_u);
Kfs_EF = Kfsbar(BC.fixed_ufdot, BC.free_u);
Kfs_FE = Kfsbar(BC.free_ufdot, BC.fixed_u);
Kfs_FF = Kfsbar(BC.free_ufdot, BC.free_u);

Kff_EE = Kffbar(BC.fixed_ufdot, BC.fixed_ufdot);
Kff_EF = Kffbar(BC.fixed_ufdot, BC.free_ufdot);
Kff_FE = Kffbar(BC.free_ufdot, BC.fixed_ufdot);
Kff_FF = Kffbar(BC.free_ufdot, BC.free_ufdot);

Kfp_EE = Kfpbar(BC.fixed_ufdot, BC.fixed_p);
Kfp_EF = Kfpbar(BC.fixed_ufdot, BC.free_p);
Kfp_FE = Kfpbar(BC.free_ufdot, BC.fixed_p);
Kfp_FF = Kfpbar(BC.free_ufdot, BC.free_p);

Kps_EE = Kpsbar(BC.fixed_p, BC.fixed_u);
Kps_EF = Kpsbar(BC.fixed_p, BC.free_u);
Kps_FE = Kpsbar(BC.free_p, BC.fixed_u);
Kps_FF = Kpsbar(BC.free_p, BC.free_u);

Kpf_EE = Kpfbar(BC.fixed_p, BC.fixed_ufdot);
Kpf_EF = Kpfbar(BC.fixed_p, BC.free_ufdot);
Kpf_FE = Kpfbar(BC.free_p, BC.fixed_ufdot);
Kpf_FF = Kpfbar(BC.free_p, BC.free_ufdot);

Kpp_EE = Kppbar(BC.fixed_p, BC.fixed_p);
Kpp_EF = Kppbar(BC.fixed_p, BC.free_p);
Kpp_FE = Kppbar(BC.free_p, BC.fixed_p);
Kpp_FF = Kppbar(BC.free_p, BC.free_p);

% matrices reassemble
KEE = [Kss_EE, Ksf_EE, Ksp_EE;
    Kfs_EE, Kff_EE, Kfp_EE;
    Kps_EE, Kpf_EE, Kpp_EE];
KEF = [Kss_EF, Ksf_EF, Ksp_EF;
    Kfs_EF, Kff_EF, Kfp_EF;
    Kps_EF, Kpf_EF, Kpp_EF];
KFE = [Kss_FE, Ksf_FE, Ksp_FE;
    Kfs_FE, Kff_FE, Kfp_FE;
    Kps_FE, Kpf_FE, Kpp_FE];
KFF = [Kss_FF, Ksf_FF, Ksp_FF;
    Kfs_FF, Kff_FF, Kfp_FF;
    Kps_FF, Kpf_FF, Kpp_FF];

% at first step: compute solid acceleration and pressure gradient
if Control.step == 1
    % find u2dot_old, uf2dot_old, pdot_old
    uf2dot_old = Mff\(ff-Kff*ufdot_old + Cfs*udot_old+Kfp*p_old);
    u2dot_old = Mss\(fu-Kss*u_old+Ksp*p_old-Msf*uf2dot_old);
    pdot_old = Kpp\(fp-Kps*udot_old - Kpf*ufdot_old);
end

% auxiliar terms for external forces vector
fubar = fu + Mss*(u_old/beta/dt^2 + udot_old/beta/dt + u2dot_old*(1/2/beta-1)) + ...
    Msf*(ufdot_old/theta/dt - uf2dot_old*(1-1/theta));

ffbar = ff + Mff*(ufdot_old/theta/dt - uf2dot_old*(1-1/theta)) -...
    Cfs*(u_old*gamma/beta/dt + udot_old*(gamma/beta-1) + u2dot_old*dt*(gamma/2/beta-1));

fpbar = fp + Kps*(u_old*gamma/beta/dt + udot_old*(gamma/beta-1) + u2dot_old*dt*(gamma/2/beta-1)) + ...
    Kpp*(p_old*(1/alpha/dt) - pdot_old*(1-1/alpha));

fuF = fubar(BC.free_u);
ffF = ffbar(BC.free_ufdot);
fpF = fpbar(BC.free_p);

fuE = fubar(BC.fixed_u);
ffE = ffbar(BC.fixed_ufdot);
fpE = fpbar(BC.fixed_p);

uE = BC.fixed_u_value(Control.t);
vE = BC.fixed_ufdot_value(Control.t);
pE = BC.fixed_p_value(Control.t);

dE = [uE; vE; pE];
fE = [fuE; ffE; fpE];
fF = [fuF; ffF; fpF];

%% Solve linear system
% solve for displacement and pressure
dF = MatrixInvert(KFF, fF - KFE *dE, Control.parallel);

% solve for reactions
rE = KEE*dE + KEF*dF - fE;

uF = dF(1:length(BC.free_u),1);
vF = dF(length(BC.free_u) + 1 : length(BC.free_u) + length(BC.free_ufdot),1);
pF = dF(length(BC.free_u) + length(BC.free_ufdot) + 1 : end,1);

%% Store u/p/uf
% force reactions
fuE = rE(1:length(BC.fixed_u),1);
% force reactions
ffE = rE(length(BC.fixed_u) + 1 : length(BC.fixed_u) + length(BC.fixed_ufdot),1);
% flux reactions
fpE = rE(length(BC.fixed_u) + length(BC.fixed_ufdot) + 1 : end,1);

u(BC.fixed_u, 1) = uE;
u(BC.free_u, 1) = uF;
ufdot(BC.fixed_ufdot, 1) = vE;
ufdot(BC.free_ufdot, 1) = vF;
p(BC.fixed_p, 1) = pE;
p(BC.free_p, 1) = pF;

%% Velocity and acceleration
udot = (u - u_old)*gamma/beta/dt - udot_old * (gamma/beta -1) - u2dot_old * dt * (gamma/2/beta-1);
u2dot = (u - u_old)/beta/dt^2 - udot_old/beta/dt - u2dot_old * (1/2/beta -1);
uf2dot = (ufdot - ufdot_old)/theta/dt + (1-1/theta)* uf2dot_old;
pdot = (p - p_old)/(lambda*dt) - (1/lambda - 1) * pdot_old;

%% Store variables
Solution.u = u;
Solution.udot = udot;
Solution.u2dot = u2dot;
Solution.p = p;
Solution.pdot = pdot;
Solution.ufdot = ufdot;
Solution.uf2dot = uf2dot;
Solution.fuE = fuE;
Solution.fpE = fpE;
Solution.ffE = ffE;

end