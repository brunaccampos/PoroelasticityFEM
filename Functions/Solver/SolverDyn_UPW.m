function [Solution] = SolverDyn_UPW(Mss, Msf, Kss, Ksp, Mfs, Mff, Kfp, Kff, Cps, Kpf, Cpp, fu, fp, ff, BC, Control, Iteration)
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
%   w: relative fluid velocity
%   wdot: relative fluid acceleration
% ------------------------------------------------------------------------

%% Iteration data
u_old = Iteration.u_old;
udot_old = Iteration.udot_old;
u2dot_old = Iteration.u2dot_old;
p_old = Iteration.p_old;
pdot_old = Iteration.pdot_old;
w_old = Iteration.w_old;
wdot_old = Iteration.wdot_old;

% current time step
dt = Control.dtc;

% time integration parameters
beta = Control.beta;
gamma = Control.gamma;
theta = Control.theta;
alpha = theta;

%% Matrix partitioning
% matrices time discretization
Kssbar = Kss + Mss./(beta*dt^2);
Ksfbar = Msf./(theta*dt);
Kspbar = -Ksp;

Kfsbar = Mfs./(beta*dt^2);
Kffbar = Mff./(theta*dt) + Kff;
Kfpbar = -Kfp;

Kpsbar = Cps*gamma/beta/dt;
Kpfbar = Kpf;
% Kppbar = zeros(length(fp), length(fp));
Kppbar = Cpp./(alpha*dt);

% matrix partitioning
Kss_EE = Kssbar(BC.fixed_u, BC.fixed_u);
Kss_EF = Kssbar(BC.fixed_u, BC.free_u);
Kss_FE = Kssbar(BC.free_u, BC.fixed_u);
Kss_FF = Kssbar(BC.free_u, BC.free_u);

Ksf_EE = Ksfbar(BC.fixed_u, BC.fixed_w);
Ksf_EF = Ksfbar(BC.fixed_u, BC.free_w);
Ksf_FE = Ksfbar(BC.free_u, BC.fixed_w);
Ksf_FF = Ksfbar(BC.free_u, BC.free_w);

Ksp_EE = Kspbar(BC.fixed_u, BC.fixed_p);
Ksp_EF = Kspbar(BC.fixed_u, BC.free_p);
Ksp_FE = Kspbar(BC.free_u, BC.fixed_p);
Ksp_FF = Kspbar(BC.free_u, BC.free_p);

Kfs_EE = Kfsbar(BC.fixed_w, BC.fixed_u);
Kfs_EF = Kfsbar(BC.fixed_w, BC.free_u);
Kfs_FE = Kfsbar(BC.free_w, BC.fixed_u);
Kfs_FF = Kfsbar(BC.free_w, BC.free_u);

Kff_EE = Kffbar(BC.fixed_w, BC.fixed_w);
Kff_EF = Kffbar(BC.fixed_w, BC.free_w);
Kff_FE = Kffbar(BC.free_w, BC.fixed_w);
Kff_FF = Kffbar(BC.free_w, BC.free_w);

Kfp_EE = Kfpbar(BC.fixed_w, BC.fixed_p);
Kfp_EF = Kfpbar(BC.fixed_w, BC.free_p);
Kfp_FE = Kfpbar(BC.free_w, BC.fixed_p);
Kfp_FF = Kfpbar(BC.free_w, BC.free_p);

Kps_EE = Kpsbar(BC.fixed_p, BC.fixed_u);
Kps_EF = Kpsbar(BC.fixed_p, BC.free_u);
Kps_FE = Kpsbar(BC.free_p, BC.fixed_u);
Kps_FF = Kpsbar(BC.free_p, BC.free_u);

Kpf_EE = Kpfbar(BC.fixed_p, BC.fixed_w);
Kpf_EF = Kpfbar(BC.fixed_p, BC.free_w);
Kpf_FE = Kpfbar(BC.free_p, BC.fixed_w);
Kpf_FF = Kpfbar(BC.free_p, BC.free_w);

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
    % [Kss 0; 0 Kff]*[u_old w_old]
    aux1 = [Kss, sparse(length(Kss), length(Kss));
        sparse(length(Kss), length(Kss)), Kff] * [u_old; w_old];
    % zeros for p 
    aux2 = [aux1; zeros(length(fp),1)];
    % find u2dot_old, wdot_old, p_old
    aux3 = [Mss, Msf, -Ksp;
        Mfs, Mff, -Kfp;
        Cps, Kpf, sparse(length(fp), length(fp))];
    aux4 = aux3 \ ([fu; ff; fp] - aux2);
    % split u2dot_old, wdot_old, p_old
    u2dot_old = aux4(1:length(u2dot_old));
    wdot_old = aux4(length(u2dot_old)+1:length(u2dot_old)+length(wdot_old));
    p_old = aux4(length(u2dot_old)+length(wdot_old)+1:end);
end

% auxiliar terms for external forces vector
fubar = fu + Mss*(u_old/beta/dt^2 + udot_old/beta/dt + (1/2/beta-1)*u2dot_old) +...
    Msf*(w_old/theta/dt - (1-1/theta)*wdot_old);

ffbar = ff + Mfs*(u_old/beta/dt^2 + udot_old/beta/dt + (1/2/beta-1)*u2dot_old) +...
    Mff*(w_old/theta/dt - (1-1/theta)*wdot_old);

fpbar = fp + Cps*(u_old*gamma/beta/dt + (gamma/beta-1)*udot_old + dt*(gamma/2/beta-1)*u2dot_old) + ...
    Cpp*(p_old/alpha/dt - (1-1/alpha)*pdot_old);

fuF = fubar(BC.free_u);
ffF = ffbar(BC.free_w);
fpF = fpbar(BC.free_p);

fuE = fubar(BC.fixed_u);
ffE = ffbar(BC.fixed_w);
fpE = fpbar(BC.fixed_p);

uE = BC.fixed_u_value(Control.t);
wE = BC.fixed_w_value(Control.t);
pE = BC.fixed_p_value(Control.t);

dE = [uE; wE; pE];
fE = [fuE; ffE; fpE];
fF = [fuF; ffF; fpF];

%% Solve linear system
% solve for displacement and pressure
dF = MatrixInvert(KFF, fF - KFE *dE, Control.parallel);

% solve for reactions
rE = KEE*dE + KEF*dF - fE;

uF = dF(1:length(BC.free_u),1);
wF = dF(length(BC.free_u)+1 : length(BC.free_u) + length(BC.free_w),1);
pF = dF(length(BC.free_u) + length(BC.free_w) + 1 : end,1);

%% Store u/p/uf
% force reactions
fuE = rE(1:length(BC.fixed_u),1);
% flux reactions
ffE = rE(length(BC.fixed_u) + 1: length(BC.fixed_u) + length(BC.fixed_w));
% force reactions
fpE = rE(length(BC.fixed_u) + length(BC.fixed_w) + 1: end,1);

u(BC.fixed_u, 1) = uE;
u(BC.free_u, 1) = uF;
w(BC.fixed_w, 1) = wE;
w(BC.free_w, 1) = wF;
p(BC.fixed_p, 1) = pE;
p(BC.free_p, 1) = pF;

%% Velocity and acceleration
udot = (u - u_old)*gamma/(beta*dt) - udot_old * (gamma/beta -1) - u2dot_old * dt * (gamma/(2*beta)-1);
u2dot = (u - u_old)/(beta*dt^2) - udot_old/(beta*dt) - u2dot_old * (1/(2*beta) -1);
pdot = (p - p_old)/(theta*dt) - (1/theta - 1) * pdot_old;
wdot = (w - w_old)/(theta*dt) + (1-1/theta)* wdot_old;

%% Store variables
Solution.u = u;
Solution.udot = udot;
Solution.u2dot = u2dot;
Solution.p = p;
Solution.pdot = pdot;
Solution.w = w;
Solution.wdot = wdot;
Solution.fuE = fuE;
Solution.fpE = fpE;
Solution.ffE = ffE;

end