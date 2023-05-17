function [Solution] = SolverDynamic_BiotUUP(Kss, Ksp, Mss, Csf, Css, Kpf, Kps, Kpp, Kfp, Mff, Cff, Cfs, fu, ff, BC, Control, Iteration)
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
if ~isempty(Iteration)
    u_old = Iteration.u_old;
    udot_old = Iteration.udot_old;
    u2dot_old = Iteration.u2dot_old;
    p_old = Iteration.p_old;
    pdot_old = Iteration.pdot_old;
    uf_old = Iteration.uf_old;
    ufdot_old = Iteration.ufdot_old;
    uf2dot_old = Iteration.uf2dot_old;
else
    u_old = zeros(length(Kss),1);
    udot_old = zeros(length(Kss),1);
    u2dot_old = zeros(length(Kss),1);
    p_old = zeros(length(Kpp),1);
    pdot_old = zeros(length(Kpp),1);
    uf_old = zeros(length(Kff),1);
    ufdot_old = zeros(length(Kff),1);
    uf2dot_old = zeros(length(Kff),1);
end

% time step
dt = Control.dt;

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
Ksfbar = -Csf*(xi/lambda/dt);

Kpsbar = Kps;
Kppbar = Kpp;
Kpfbar = Kpf;

Kfsbar = -Cfs*(gamma/beta/dt);
Kfpbar = -Kfp;
Kffbar = Mff*(1/lambda/dt^2) + Cff*(xi/lambda/dt);

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

% matrices for unknown DOFs
MssFF = Mss(BC.free_u, BC.free_u);
KssFF = Kss(BC.free_u, BC.free_u);
KspFF = Ksp(BC.free_u, BC.free_p);
CsfFF = Csf(BC.free_u, BC.free_u);
CssFF = Css(BC.free_u, BC.free_u);

MffFF = Mff(BC.free_u, BC.free_u);
KfpFF = Kfp(BC.free_u, BC.free_p);
CffFF = Cff(BC.free_u, BC.free_u);
CfsFF = Cfs(BC.free_u, BC.free_u);

% at first step: compute solid acceleration and pressure gradient
if Control.step == 1
    u2dot_old(BC.free_u) = MssFF\(fu(BC.free_u) - KssFF*u_old(BC.free_u) + ...
        KspFF*p_old(BC.free_p) + CsfFF*ufdot_old(BC.free_u) - CssFF*udot_old(BC.free_u));
    uf2dot_old(BC.free_u) = MffFF\(ff(BC.free_u) + KfpFF*p_old(BC.free_p) - ...
        CffFF*ufdot_old(BC.free_u) + CfsFF*udot_old(BC.free_u));
end

% auxiliar terms for external forces vector
fuF = fu(BC.free_u) + MssFF * (1/(beta*dt^2) * u_old(BC.free_u) + 1/(beta*dt) * udot_old(BC.free_u) + ...
    (1/(2*beta) -1) * u2dot_old(BC.free_u)) + CsfFF * (-xi/(lambda*dt) * uf_old(BC.free_u) - ...
    (xi/lambda-1) * ufdot_old(BC.free_u) - dt*(xi/(2*lambda)-1) * uf2dot_old(BC.free_u)) + ...
    CssFF * (gamma/(beta*dt) * u_old(BC.free_u) +(gamma/beta-1) * udot_old(BC.free_u) + dt*(gamma/(2*lambda)-1) * u2dot_old(BC.free_u));

ffF = ff(BC.free_u) + MffFF * (1/(lambda*dt^2) * uf_old(BC.free_u) + 1/(lambda*dt)*ufdot_old(BC.free_u) + ...
    (1/(2*lambda)-1) * uf2dot_old(BC.free_u)) + CffFF * (xi/(lambda*dt) * uf_old(BC.free_u) + ...
    (xi/lambda-1) * ufdot_old(BC.free_u) + dt*(xi/(2*lambda)-1) * uf2dot_old(BC.free_u)) + ...
    CfsFF * (-gamma/(beta*dt) * u_old(BC.free_u) - (gamma/beta-1) * udot_old(BC.free_u) - dt*(gamma/(2*beta)-1) * u2dot_old(BC.free_u));


fuE = fu(BC.fixed_u);
fpE = zeros(length(BC.fixed_p),1);
ffE = ff(BC.fixed_u);

uE = BC.fixed_u_value;
pE = zeros(length(BC.fixed_p),1);
ufE = BC.fixed_u_value;

dE = [uE; pE; ufE];
fE = [fuE; fpE; ffE];
fF = [fuF; zeros(length(BC.free_p),1); ffF];

%% Solve linear system
% solve for displacement and pressure
dF = KFF\(fF - KFE *dE);

% solve for reactions
rE = KEE*dE + KEF*dF - fE;

uF = dF(1:length(BC.free_u),1);
pF = dF(length(BC.free_u)+1 : length(BC.free_u) + length(BC.free_p),1);
ufF = dF(length(BC.free_u) + length(BC.free_p) + 1 : end,1);

%% Store u/p/uf
% force reactions
fuE = rE(1:length(BC.fixed_u),1);
% flux reactions
fpE = rE(length(BC.fixed_u)+1 : length(BC.fixed_u) + length(BC.fixed_p),1);
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