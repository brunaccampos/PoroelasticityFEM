function [Solution] = SolverTr_UPN(Kuu, Kup, Kpp, Kpu, S, Kpn, Knn, Knu, Knp, Kun, fu, fp, fn, BC, Control, Iteration)
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
% version 6: changing time discretization to beta method
% ------------------------------------------------------------------------

%% Iteration data
if ~isempty(Iteration)
    u_old = Iteration.u_old;
    p_old = Iteration.p_old;
    n_old = Iteration.n_old;
    fu_old = Iteration.fu_old;
    fp_old = Iteration.fp_old;
    fn_old = Iteration.fn_old;
else
    u_old = zeros(length(Kuu),1);
    p_old = zeros(length(Kpp),1);
    n_old = zeros(length(Knn),1);
    fu_old = zeros(length(Kuu),1);
    fp_old = zeros(length(Kpp),1);
    fn_old = zeros(length(Knn),1);
end
% time step
dt = Control.dt;

% beta parameter
beta = Control.beta;

%% Matrix partitioning
% matrices time discretization
Kuubar = beta * dt * Kuu;
Kupbar = -beta * dt * Kup;
Kunbar = Kun;

Kpubar = Kpu;
Kppbar = beta * dt * Kpp + S;
Kpnbar = Kpn;

Knubar = Knu;
Knpbar = beta * dt * Knp;
Knnbar = Knn;


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

% auxiliar terms for external forces vector
fubar = -dt * (1-beta) * Kuu * u_old + dt * (1-beta) * Kup * p_old + dt * (1-beta) * fu_old + dt * beta * fu;
fpbar = Kpu * u_old + (S - dt * (1-beta) * Kpp) * p_old + Kpn * n_old + dt * (1-beta) * fp_old + dt * beta * fp;
fnbar = Knu * u_old - dt * (1-beta) * Knp * p_old + Knn * n_old + dt * (1-beta) * fn_old + dt * beta * fn;


% partitioning vectors
fuF = fubar(BC.free_u);
fpF = fpbar(BC.free_p);
fnF = fnbar(BC.free_n);

fuE = fubar(BC.fixed_u);
fpE = fpbar(BC.fixed_p);
fnE = fnbar(BC.fixed_n);

uE = BC.fixed_u_value(Control.t);
pE = BC.fixed_p_value(Control.t);
nE = BC.fixed_n_value;

dE = [uE; pE; nE];
fF = [fuF; fpF; fnF];
fE = [fuE; fpE; fnE];

%% Solve linear system
% solve for displacement and pressure
dF = MatrixInvert(KFF, fF - KFE *dE, Control.parallel);

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
fnE = rE(length(BC.fixed_u) + length(BC.fixed_p) + 1:end,1);

u(BC.fixed_u, 1) = uE;
u(BC.free_u, 1) = uF;
p(BC.fixed_p, 1) = pE;
p(BC.free_p, 1) = pF;
n(BC.fixed_n, 1) = nE;
n(BC.free_n, 1) = nF;

%% Gradients
udot = (u - u_old)./(beta*dt);
pdot = (p - p_old)./(beta*dt);
ndot = (n - n_old)./(beta*dt);

%% Store variables
Solution.u = u;
Solution.udot = udot;
Solution.p = p;
Solution.pdot = pdot;
Solution.n = n;
Solution.ndot = ndot;
Solution.fE = fE;
Solution.qE = qE;
Solution.fnE = fnE;

end