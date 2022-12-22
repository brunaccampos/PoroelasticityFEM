function [Solution] = SolverDynamic_Spanos(Muu, Mpu, Mnu, Kuu, Kup, Kpp, Kpu, S, Kpn, Knn, Knu, Knp, Kun, fu, fp, fn, BC, Control, Iteration)
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
% ------------------------------------------------------------------------
% version 5: correct Spanos formulation
% ------------------------------------------------------------------------

%% Iteration data
if ~isempty(Iteration)
    u_old = Iteration.u_old;
    udot_old = Iteration.udot_old;
    u2dot_old = Iteration.u2dot_old;
    p_old = Iteration.p_old;
    pdot_old = Iteration.pdot_old;
    n_old = Iteration.n_old;
    ndot_old = Iteration.ndot_old;
else
    u_old = zeros(length(Kuu),1);
    udot_old = zeros(length(Kuu),1);
    u2dot_old = zeros(length(Kuu),1);
    p_old = zeros(length(Kpp),1);
    pdot_old = zeros(length(Kpp),1);
    n_old = zeros(length(Knn),1);
    ndot_old = zeros(length(Knn),1);
end

% time step
dt = Control.dt;

% time integration parameters
beta = Control.beta;
gamma = Control.gamma;
theta = Control.theta;
lambda = Control.lambda;

%% Matrix partitioning
% matrices time discretization
Kuubar = Kuu + Muu./(beta*dt^2);
Kupbar = -Kup;
Kunbar = Kun;

Kpubar = Kpu *gamma/(beta*dt) - Mpu/(beta*dt^2);
Kppbar = Kpp + S/(theta*dt);
Kpnbar = Kpn./(lambda*dt);

Knubar = Knu *gamma/(beta*dt) - Mnu./(beta*dt^2);
Knpbar = Knp;
Knnbar = Knn./(lambda*dt);

% matrix partitioning
Kuu_EE = Kuubar(BC.fixed_u, BC.fixed_u);
Kuu_EF = Kuubar(BC.fixed_u, BC.free_u);
Kuu_FE = Kuubar(BC.free_u, BC.fixed_u);
Kuu_FF = Kuubar(BC.free_u, BC.free_u);

Kup_EE = Kupbar(BC.fixed_u, BC.fixed_p);
Kup_EF = Kupbar(BC.fixed_u, BC.free_p);
Kup_FE = Kupbar(BC.free_u, BC.fixed_p);
Kup_FF = Kupbar(BC.free_u, BC.free_p);

Kpu_EE = Kpubar(BC.fixed_p, BC.fixed_u);
Kpu_EF = Kpubar(BC.fixed_p, BC.free_u);
Kpu_FE = Kpubar(BC.free_p, BC.fixed_u);
Kpu_FF = Kpubar(BC.free_p, BC.free_u);

Kpp_EE = Kppbar(BC.fixed_p, BC.fixed_p);
Kpp_EF = Kppbar(BC.fixed_p, BC.free_p);
Kpp_FE = Kppbar(BC.free_p, BC.fixed_p);
Kpp_FF = Kppbar(BC.free_p, BC.free_p);

Kun_EE = Kunbar(BC.fixed_u, BC.fixed_n);
Kun_EF = Kunbar(BC.fixed_u, BC.free_n);
Kun_FE = Kunbar(BC.free_u, BC.fixed_n);
Kun_FF = Kunbar(BC.free_u, BC.free_n);

Kpn_EE = Kpnbar(BC.fixed_p, BC.fixed_n);
Kpn_EF = Kpnbar(BC.fixed_p, BC.free_n);
Kpn_FE = Kpnbar(BC.free_p, BC.fixed_n);
Kpn_FF = Kpnbar(BC.free_p, BC.free_n);

Knu_EE = Knubar(BC.fixed_n, BC.fixed_u);
Knu_EF = Knubar(BC.fixed_n, BC.free_u);
Knu_FE = Knubar(BC.free_n, BC.fixed_u);
Knu_FF = Knubar(BC.free_n, BC.free_u);

Knp_EE = Knpbar(BC.fixed_n, BC.fixed_p);
Knp_EF = Knpbar(BC.fixed_n, BC.free_p);
Knp_FE = Knpbar(BC.free_n, BC.fixed_p);
Knp_FF = Knpbar(BC.free_n, BC.free_p);

Knn_EE = Knnbar(BC.fixed_n, BC.fixed_n);
Knn_EF = Knnbar(BC.fixed_n, BC.free_n);
Knn_FE = Knnbar(BC.free_n, BC.fixed_n);
Knn_FF = Knnbar(BC.free_n, BC.free_n);

% matrices reassemble
KEE = [Kuu_EE, Kup_EE, Kun_EE;
    Kpu_EE, Kpp_EE, Kpn_EE;
    Knu_EE, Knp_EE, Knn_EE];
KEF = [Kuu_EF, Kup_EF, Kun_EF;
    Kpu_EF, Kpp_EF, Kpn_EF;
    Knu_EF, Knp_EF, Knn_EF];
KFE = [Kuu_FE, Kup_FE, Kun_FE;
    Kpu_FE, Kpp_FE, Kpn_FE;
    Knu_FE, Knp_FE, Knn_FE];
KFF = [Kuu_FF, Kup_FF, Kun_FF;
    Kpu_FF, Kpp_FF, Kpn_FF;
    Knu_FF, Knp_FF, Knn_FF];

% matrices for unknown DOFs
MuuFF = Muu(BC.free_u, BC.free_u);
KuuFF = Kuu(BC.free_u, BC.free_u);
KupFF = Kup(BC.free_u, BC.free_p);

MpuFF = Mpu(BC.free_p, BC.free_u);
KpuFF = Kpu(BC.free_p, BC.free_u);
KppFF = Kpp(BC.free_p, BC.free_p);
SFF = S(BC.free_p, BC.free_p);
KpnFF = Kpn(BC.free_p, BC.free_n);

% at first step: compute solid acceleration and pressure gradient
if Control.step == 1
    u2dot_old(BC.free_u) = MuuFF\(-KuuFF*u_old(BC.free_u) + KupFF*p_old(BC.free_p) + fu(BC.free_u));
    pdot_old(BC.free_p) = SFF\(MpuFF*u2dot_old(BC.free_u) - KpuFF*udot_old(BC.free_u) - KppFF*p_old(BC.free_p) - KpnFF*ndot_old(BC.free_n) + fp(BC.free_p));
end

% auxiliar terms for external forces vector
fubar = fu + Muu * (u_old./(beta*dt^2) + udot_old./(beta*dt) + (1/(2*beta)-1) * u2dot_old);
fpbar = fp - Mpu * (u_old./(beta*dt^2) + udot_old./(beta*dt) + (1/(2*beta)-1) * u2dot_old) + ...
    Kpu * (u_old * gamma/(beta*dt) + (gamma/beta -1) * udot_old + dt*(gamma/(2*beta)-1) * u2dot_old) + ...
    S * (p_old/(theta*dt) + (1/theta-1) * pdot_old) + Kpn * (n_old/(lambda*dt) + (1/lambda-1) * ndot_old);
fnbar = fn - Mnu * (u_old./(beta*dt^2) + udot_old ./ (beta*dt) + (1/(2*beta)-1) * u2dot_old) + ...
    Knu * (u_old *gamma/(beta*dt) + (gamma/beta-1) * udot_old + dt*(gamma/(2*beta)-1) * u2dot_old) + ...
    Knn * (n_old./(lambda*dt) + (1/lambda-1) * ndot_old);

% partitioning vectors
fuF = fubar(BC.free_u);
fpF = fpbar(BC.free_p);
fnF = fnbar(BC.free_n);

fuE = fubar(BC.fixed_u);
fpE = fpbar(BC.fixed_p);
fnE = fnbar(BC.fixed_n);

uE = BC.fixed_u_value;
pE = BC.fixed_p_value;
nE = BC.fixed_n_value;

dE = [uE; pE; nE];
fF = [fuF; fpF; fnF];
fE = [fuE; fpE; fnE];

%% Solve linear system
% solve for displacement and pressure
dF = KFF\(fF - KFE *dE);

% solve for reactions
rE = KEE*dE + KEF*dF - fE;

%% Store u/p/n
uF = dF(1:length(BC.free_u),1);
pF = dF(length(BC.free_u)+1 : length(BC.free_u) + length(BC.free_p),1);
nF = dF(length(BC.free_u) + length(BC.free_p) + 1 : end,1);

% force reactions
fE = rE(1:length(BC.fixed_u),1);
% flux reactions
qE = rE(length(BC.fixed_u)+1 : length(BC.fixed_u) + length(BC.fixed_p),1);

u(BC.fixed_u, 1) = uE;
u(BC.free_u, 1) = uF;
p(BC.fixed_p, 1) = pE;
p(BC.free_p, 1) = pF;
n(BC.fixed_n, 1) = nE;
n(BC.free_n, 1) = nF;

%% Velocity and acceleration
udot = (u - u_old)*gamma/(beta*dt) - udot_old * (gamma/beta -1) - u2dot_old * dt * (gamma/(2*beta)-1);
u2dot = (u - u_old)/(beta*dt^2) - udot_old/(beta*dt) - u2dot_old * (1/(2*beta) -1);
pdot = (p - p_old)/(theta*dt) - (1/theta - 1) * pdot_old;
ndot = (n - n_old)/(lambda*dt) - (1/lambda - 1) * ndot_old;

%% Store variables
Solution.u = u;
Solution.udot = udot;
Solution.u2dot = u2dot;
Solution.p = p;
Solution.pdot = pdot;
Solution.n = n;
Solution.ndot = ndot;
Solution.fE = fE;
Solution.qE = qE;

end