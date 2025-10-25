% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function [Solution] = SolverTr_UPU(Kss, Ksp, Csf, Css, Kpf, Kps, Kpp, Kfp, Cff, Cfs, fu, fp, ff, BC, Control, Iteration)
% Solve linear system for transient case

%% Iteration data
u_old = Iteration.u_old;
p_old = Iteration.p_old;
uf_old = Iteration.uf_old;
fs_old = Iteration.fu_old;
fp_old = Iteration.fp_old;
ff_old = Iteration.ff_old;

% current time step
dt = Control.dtc;

% beta parameter
beta = Control.beta;

%% Matrix partitioning
% matrices time discretization
Kssbar = Css + dt*beta*Kss;
Kspbar = -dt*beta*Ksp;
Ksfbar = -Csf;

Kpsbar = dt*beta*Kps;
Kppbar = dt*beta*Kpp;
Kpfbar = dt*beta*Kpf;

Kfsbar = -Cfs;
Kfpbar = -dt*beta*Kfp;
Kffbar = Cff;

% matrix partitioning
Kss_EE = Kssbar(BC.fixed_u, BC.fixed_u);
Kss_EF = Kssbar(BC.fixed_u, BC.free_u);
Kss_FE = Kssbar(BC.free_u, BC.fixed_u);
Kss_FF = Kssbar(BC.free_u, BC.free_u);

Ksp_EE = Kspbar(BC.fixed_u, BC.fixed_p);
Ksp_EF = Kspbar(BC.fixed_u, BC.free_p);
Ksp_FE = Kspbar(BC.free_u, BC.fixed_p);
Ksp_FF = Kspbar(BC.free_u, BC.free_p);

Ksf_EE = Ksfbar(BC.fixed_u, BC.fixed_u);
Ksf_EF = Ksfbar(BC.fixed_u, BC.free_u);
Ksf_FE = Ksfbar(BC.free_u, BC.fixed_u);
Ksf_FF = Ksfbar(BC.free_u, BC.free_u);

Kps_EE = Kpsbar(BC.fixed_p, BC.fixed_u);
Kps_EF = Kpsbar(BC.fixed_p, BC.free_u);
Kps_FE = Kpsbar(BC.free_p, BC.fixed_u);
Kps_FF = Kpsbar(BC.free_p, BC.free_u);

Kpp_EE = Kppbar(BC.fixed_p, BC.fixed_p);
Kpp_EF = Kppbar(BC.fixed_p, BC.free_p);
Kpp_FE = Kppbar(BC.free_p, BC.fixed_p);
Kpp_FF = Kppbar(BC.free_p, BC.free_p);

Kpf_EE = Kpfbar(BC.fixed_p, BC.fixed_u);
Kpf_EF = Kpfbar(BC.fixed_p, BC.free_u);
Kpf_FE = Kpfbar(BC.free_p, BC.fixed_u);
Kpf_FF = Kpfbar(BC.free_p, BC.free_u);

Kfs_EE = Kfsbar(BC.fixed_u, BC.fixed_u);
Kfs_EF = Kfsbar(BC.fixed_u, BC.free_u);
Kfs_FE = Kfsbar(BC.free_u, BC.fixed_u);
Kfs_FF = Kfsbar(BC.free_u, BC.free_u);

Kfp_EE = Kfpbar(BC.fixed_u, BC.fixed_p);
Kfp_EF = Kfpbar(BC.fixed_u, BC.free_p);
Kfp_FE = Kfpbar(BC.free_u, BC.fixed_p);
Kfp_FF = Kfpbar(BC.free_u, BC.free_p);

Kff_EE = Kffbar(BC.fixed_u, BC.fixed_u);
Kff_EF = Kffbar(BC.fixed_u, BC.free_u);
Kff_FE = Kffbar(BC.free_u, BC.fixed_u);
Kff_FF = Kffbar(BC.free_u, BC.free_u);

% matrices reassemble
KEE = [Kss_EE, Ksp_EE, Ksf_EE;
    Kps_EE, Kpp_EE, Kpf_EE;
    Kfs_EE, Kfp_EE, Kff_EE];
KEF = [Kss_EF, Ksp_EF, Ksf_EF;
    Kps_EF, Kpp_EF, Kpf_EF;
    Kfs_EF, Kfp_EF, Kff_EF];
KFE = [Kss_FE, Ksp_FE, Ksf_FE;
    Kps_FE, Kpp_FE, Kpf_FE;
    Kfs_FE, Kfp_FE, Kff_FE];
KFF = [Kss_FF, Ksp_FF, Ksf_FF;
    Kps_FF, Kpp_FF, Kpf_FF;
    Kfs_FF, Kfp_FF, Kff_FF];

% auxiliar terms for external forces vector
fsbar = (Css - dt*(1-beta)*Kss)*u_old + dt*(1-beta)*Ksp*p_old - Csf*uf_old +...
    dt*(1-beta)*fs_old + dt*beta*fu;

fpbar = -dt*(1-beta)*Kps*u_old - dt*(1-beta)*Kpp*p_old - dt*(1-beta)*Kpf*uf_old +...
    dt*(1-beta)*fp_old + dt*beta*fp;

ffbar = -Cfs*u_old + dt*(1-beta)*Kfp*p_old + Cff*uf_old + dt*(1-beta)*ff_old +...
    dt*beta*ff;

fsF = fsbar(BC.free_u);
fpF = fpbar(BC.free_p);
ffF = ffbar(BC.free_uf);


fuE = fsbar(BC.fixed_u);
fpE = fpbar(BC.fixed_p);
ffE = ffbar(BC.fixed_uf);

uE = BC.fixed_u_value(Control.t);
pE = BC.fixed_p_value(Control.t);
ufE = BC.fixed_uf_value(Control.t);

dE = [uE; pE; ufE];
fE = [fuE; fpE; ffE];
fF = [fsF; fpF; ffF];

%% Solve linear system
% solve for displacement and pressure
dF = MatrixInvert(KFF, fF - KFE *dE, Control.parallel);

% solve for reactions
rE = KEE*dE + KEF*dF - fE;

uF = dF(1:length(BC.free_u),1);
pF = dF(length(BC.free_u)+1 : length(BC.free_u) + length(BC.free_p),1);
ufF = dF(length(BC.free_u) + length(BC.free_p) + 1 : end,1);

%% Store u/p/uf
% force reactions
fuE = rE(1:length(BC.fixed_u),1);
% flux reactions
fpE = rE(length(BC.fixed_u) + 1: length(BC.fixed_u) + length(BC.fixed_p));
% force reactions
ffE = rE(length(BC.fixed_u) + length(BC.fixed_p) + 1: end,1);

u(BC.fixed_u, 1) = uE;
u(BC.free_u, 1) = uF;
p(BC.fixed_p, 1) = pE;
p(BC.free_p, 1) = pF;
uf(BC.fixed_uf, 1) = ufE;
uf(BC.free_uf, 1) = ufF;

%% Velocity and acceleration
udot = (u - u_old)./(beta*dt);
pdot = (p - p_old)./(beta*dt);
ufdot = (uf - uf_old)./(beta*dt);

% uncomment for velocity impact problem
% udot(BC.fixed_u2) = 1;
% ufdot(BC.fixed_u2) = 1;

%% Store variables
Solution.u = u;
Solution.udot = udot;
Solution.p = p;
Solution.pdot = pdot;
Solution.uf = uf;
Solution.ufdot = ufdot;
Solution.fuE = fuE;
Solution.fpE = fpE;
Solution.ffE = ffE;

end