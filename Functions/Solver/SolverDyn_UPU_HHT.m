function [Solution] = SolverDyn_UPU_HHT(Kss, Ksp, Mss, Csf, Css, Kpf, Kps, Kpp, Kfp, Mff, Cff, Cfs, fu, fp, ff, BC, Control, Iteration)
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
%   uf: fluid displacement
%   ufdot: fluid velocity
%   uf2dot: fluid acceleration
% ------------------------------------------------------------------------

%% Iteration data
u_old = Iteration.u_old;
udot_old = Iteration.udot_old;
u2dot_old = Iteration.u2dot_old;
p_old = Iteration.p_old;
pdot_old = Iteration.pdot_old;
uf_old = Iteration.uf_old;
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

% HHT method parameter
alpha = Control.alpha;

%% Matrix partitioning
% matrices time discretization
Kssbar = (1+alpha)*Kss + Mss./(beta*dt^2) + (1+alpha)*Css*(gamma/beta/dt);
Kspbar = -(1+alpha)*Ksp;
Ksfbar = -(1+alpha)*Csf*(xi/lambda/dt);

Kpsbar = -(1+alpha)*Kps;
Kppbar = -(1+alpha)*Kpp;
Kpfbar = -(1+alpha)*Kpf;

Kfsbar = -(1+alpha)*Cfs*(gamma/beta/dt);
Kfpbar = -(1+alpha)*Kfp;
Kffbar = Mff*(1/lambda/dt^2) + (1+alpha)*Cff*(xi/lambda/dt);

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

% at first step: compute solid acceleration and pressure gradient
if Control.step == 1
    u2dot_old = Mss\(fu - Kss*u_old + Ksp*p_old + Csf*ufdot_old - Css*udot_old);
    uf2dot_old = Mff\(ff + Kfp*p_old - Cff*ufdot_old + Cfs*udot_old);
end

% auxiliar terms for external forces vector
fubar = fu + Mss*(1/(beta*dt^2) * u_old + 1/(beta*dt) * udot_old + (1/(2*beta) -1) * u2dot_old) +...
    (1+alpha)*Csf* (-xi/(lambda*dt) * uf_old - (xi/lambda-1) * ufdot_old - dt*(xi/(2*lambda)-1) * uf2dot_old) -...
    alpha*Csf * ufdot_old + (1+alpha)*Css * (gamma/(beta*dt) * u_old +(gamma/beta-1) * udot_old + ...
    dt*(gamma/(2*lambda)-1) * u2dot_old) + alpha*Css * udot_old + alpha*Kss * u_old - ...
    alpha*Ksp * p_old;

fpbar = fp - alpha*Kps * u_old - alpha*Kpp * p_old - alpha*Kpf * uf_old;

ffbar = ff + Mff * (1/(lambda*dt^2) * uf_old + 1/(lambda*dt)*ufdot_old + ...
    (1/(2*lambda)-1) * uf2dot_old) + (1+alpha)*Cff * (xi/(lambda*dt) * uf_old + ...
    (xi/lambda-1) * ufdot_old + dt*(xi/(2*lambda)-1) * uf2dot_old) + ...
    (1+alpha)*Cfs * (-gamma/(beta*dt) * u_old - (gamma/beta-1) * udot_old - ...
    dt*(gamma/(2*beta)-1) * u2dot_old) + alpha*Cff * ufdot_old - alpha*Cfs * udot_old - ...
    alpha*Kfp * p_old;

fuF = fubar(BC.free_u);
fpF = fpbar(BC.free_p);
ffF = ffbar(BC.free_uf);

fuE = fubar(BC.fixed_u);
fpE = fpbar(BC.fixed_p);
ffE = ffbar(BC.fixed_uf);

uE = BC.fixed_u_value(Control.t);
pE = BC.fixed_p_value(Control.t);
ufE = BC.fixed_uf_value(Control.t);

% uncomment for velocity impact problem (1D)
% uE(end) = uE(end)*Control.t;
% ufE(end) = ufE(end)*Control.t;

dE = [uE; pE; ufE];
fE = [fuE; fpE; ffE];
fF = [fuF; fpF; ffF];

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
% force reactions
ffE = rE(length(BC.fixed_u) + length(BC.fixed_u) + 1: end,1);

u(BC.fixed_u, 1) = uE;
u(BC.free_u, 1) = uF;
p(BC.fixed_p, 1) = pE;
p(BC.free_p, 1) = pF;
uf(BC.fixed_u, 1) = ufE;
uf(BC.free_u, 1) = ufF;

%% Velocity and acceleration
udot = (u - u_old)*gamma/(beta*dt) - udot_old * (gamma/beta -1) - u2dot_old * dt * (gamma/(2*beta)-1);
u2dot = (u - u_old)/(beta*dt^2) - udot_old/(beta*dt) - u2dot_old * (1/(2*beta) -1);
pdot = (p - p_old)/(theta*dt) - (1/theta - 1) * pdot_old;
ufdot = (uf - uf_old)*xi/(lambda*dt) - ufdot_old * (xi/lambda -1) - uf2dot_old * dt * (xi/(2*lambda)-1);
uf2dot = (uf - uf_old)/(lambda*dt^2) - ufdot_old/(lambda*dt) - uf2dot_old * (1/(2*lambda) -1);

% uncomment for velocity impact problem
% udot(BC.fixed_u2) = 1;
% ufdot(BC.fixed_u2) = 1;

%% Store variables
Solution.u = u;
Solution.udot = udot;
Solution.u2dot = u2dot;
Solution.p = p;
Solution.pdot = pdot;
Solution.uf = uf;
Solution.ufdot = ufdot;
Solution.uf2dot = uf2dot;
Solution.fuE = fuE;
Solution.ffE = ffE;

end