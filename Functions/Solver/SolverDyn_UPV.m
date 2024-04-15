function [Solution] = SolverDyn_UPV(Kss, Ksp, Mss, Csf, Css, Kpf, Kps, Kpp, Kfp, Mff, Cff, Cfs, Msf, Mfs, fu, fp, ff, BC, Control, Iteration)
% ------------------------------------------------------------------------
% Solve linear system for dynamic case
% Input parameters: coupled matrices, BC, Control, Iteration
% ------------------------------------------------------------------------
% Outputs:
%   u: solid displacement
%   udot: solid velocity
%   u2dot: solid acceleration
%   p: fluid pressure
%   pdot: fluid pressure gradient
%   ufdot: fluid velocity
%   uf2dot: fluid acceleration
% ------------------------------------------------------------------------

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
xi = gamma;

%% Matrix partitioning
% matrices time discretization
Kssbar = Kss + Mss./(beta*dt^2) + Css*(gamma/beta/dt);
Kspbar = -Ksp;
Ksfbar = Msf./(lambda*dt^2) - Csf;

Kpsbar = Kps*gamma/beta/dt;
Kppbar = Kpp/theta/dt;
Kpfbar = Kpf;

Kfsbar = Mfs./(beta*dt^2) - Cfs*(gamma/beta/dt);
Kfpbar = -Kfp;
Kffbar = Mff./(lambda*dt^2) + Cff;

% matrix partitioning
Kss_EE = Kssbar(BC.fixed_u, BC.fixed_u);
Kss_EF = Kssbar(BC.fixed_u, BC.free_u);
Kss_FE = Kssbar(BC.free_u, BC.fixed_u);
Kss_FF = Kssbar(BC.free_u, BC.free_u);

Ksp_EE = Kspbar(BC.fixed_u, BC.fixed_p);
Ksp_EF = Kspbar(BC.fixed_u, BC.free_p);
Ksp_FE = Kspbar(BC.free_u, BC.fixed_p);
Ksp_FF = Kspbar(BC.free_u, BC.free_p);

Ksf_EE = Ksfbar(BC.fixed_u, BC.fixed_ufdot);
Ksf_EF = Ksfbar(BC.fixed_u, BC.free_ufdot);
Ksf_FE = Ksfbar(BC.free_u, BC.fixed_ufdot);
Ksf_FF = Ksfbar(BC.free_u, BC.free_ufdot);

Kps_EE = Kpsbar(BC.fixed_p, BC.fixed_u);
Kps_EF = Kpsbar(BC.fixed_p, BC.free_u);
Kps_FE = Kpsbar(BC.free_p, BC.fixed_u);
Kps_FF = Kpsbar(BC.free_p, BC.free_u);

Kpp_EE = Kppbar(BC.fixed_p, BC.fixed_p);
Kpp_EF = Kppbar(BC.fixed_p, BC.free_p);
Kpp_FE = Kppbar(BC.free_p, BC.fixed_p);
Kpp_FF = Kppbar(BC.free_p, BC.free_p);

Kpf_EE = Kpfbar(BC.fixed_p, BC.fixed_ufdot);
Kpf_EF = Kpfbar(BC.fixed_p, BC.free_ufdot);
Kpf_FE = Kpfbar(BC.free_p, BC.fixed_ufdot);
Kpf_FF = Kpfbar(BC.free_p, BC.free_ufdot);

Kfs_EE = Kfsbar(BC.fixed_ufdot, BC.fixed_u);
Kfs_EF = Kfsbar(BC.fixed_ufdot, BC.free_u);
Kfs_FE = Kfsbar(BC.free_ufdot, BC.fixed_u);
Kfs_FF = Kfsbar(BC.free_ufdot, BC.free_u);

Kfp_EE = Kfpbar(BC.fixed_ufdot, BC.fixed_p);
Kfp_EF = Kfpbar(BC.fixed_ufdot, BC.free_p);
Kfp_FE = Kfpbar(BC.free_ufdot, BC.fixed_p);
Kfp_FF = Kfpbar(BC.free_ufdot, BC.free_p);

Kff_EE = Kffbar(BC.fixed_ufdot, BC.fixed_ufdot);
Kff_EF = Kffbar(BC.fixed_ufdot, BC.free_ufdot);
Kff_FE = Kffbar(BC.free_ufdot, BC.fixed_ufdot);
Kff_FF = Kffbar(BC.free_ufdot, BC.free_ufdot);

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

% at first step: compute solid acceleration and pressure gradient
if Control.step == 1
    aux = [Mss, Msf; Mfs, Mff] \([fu; ff] - [Css, -Csf; -Cfs, Cff] * [udot_old; ufdot_old]...
        - [Kss, -Ksp; sparse(length(Kss), length(Kss)), -Kfp] * [u_old; p_old]);
    u2dot_old = aux(1:length(u2dot_old));
    uf2dot_old = aux(length(u2dot_old)+1:end);
end

% auxiliar terms for external forces vector
fubar = fu + (Mss/beta/dt^2 + Css*gamma/beta*dt) * u_old + (Mss/beta/dt + Css*(gamma/beta-1)) * udot_old + ...
    (Mss*(1/2/beta-1) + Css*dt*(gamma/2/beta-1)) * u2dot_old + Msf/xi/dt * ufdot_old - Msf*(1-1/xi) * uf2dot_old;

fpbar = fp + Kps*(gamma/beta/dt*u_old + (gamma/beta-1)*udot_old + dt*(gamma/2/beta-1)*u2dot_old) + ...
    Kpp*(1/theta/dt*p_old + (1/theta-1)*pdot_old);

ffbar = ff + Mff*(1/xi/dt*ufdot_old - (1-1/xi)*uf2dot_old) + (Mfs/beta/dt^2 - Cfs*gamma/beta/dt) * u_old + ...
    (Mfs/beta/dt + Cfs*(gamma/beta-1)) * udot_old + (Mfs*(1/2/beta-1) + Cfs*dt*(gamma/2/beta-1)) * u2dot_old;

fuF = fubar(BC.free_u);
fpF = fpbar(BC.free_p);
ffF = ffbar(BC.free_uf);

fuE = fubar(BC.fixed_u);
fpE = fpbar(BC.fixed_p);
ffE = ffbar(BC.fixed_uf);

uE = BC.fixed_u_value(Control.t);
pE = BC.fixed_p_value(Control.t);
vE = BC.fixed_uf_value(Control.t);

dE = [uE; pE; vE];
fE = [fuE; fpE; ffE];
fF = [fuF; fpF; ffF];

%% Solve linear system
% solve for displacement and pressure
dF = MatrixInvert(KFF, fF - KFE *dE, Control.parallel);

% solve for reactions
rE = KEE*dE + KEF*dF - fE;

uF = dF(1:length(BC.free_u),1);
pF = dF(length(BC.free_u)+1 : length(BC.free_u) + length(BC.free_p),1);
vF = dF(length(BC.free_u) + length(BC.free_p) + 1 : end,1);

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
ufdot(BC.fixed_uf, 1) = vE;
ufdot(BC.free_uf, 1) = vF;

%% Velocity and acceleration
udot = (u - u_old)*gamma/(beta*dt) - udot_old * (gamma/beta -1) - u2dot_old * dt * (gamma/(2*beta)-1);
u2dot = (u - u_old)/(beta*dt^2) - udot_old/(beta*dt) - u2dot_old * (1/(2*beta) -1);
pdot = (p - p_old)/(theta*dt) - (1/theta - 1) * pdot_old;
uf2dot = (ufdot - ufdot_old)/(xi*dt) + (1-1/xi)* uf2dot_old;

% uncomment for velocity impact problem
% udot(BC.fixed_u2) = 1;
% ufdot(BC.fixed_u2) = 1;

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

Solution.uf = BC.initUf;

end