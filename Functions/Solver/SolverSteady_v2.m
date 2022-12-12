function [u, udot, p, pdot, n, ndot, fE, qE] = SolverSteady_v2(Kuu, Kup, Kpp, S, Kpu, Kun, Kpn, Knn, Knu, Knp, fu, fp, BC, Control, Iteration)
% Solve linear system for quasi-steady case 
% Input parameters: coupled matrices, BC, Control, Iteration
% Outputs:
%   u: solid displacement
%   udot: solid velocity
%   p: fluid pressure

%% Iteration data
u_old = Iteration.u_old;
p_old = Iteration.p_old;
n_old = Iteration.n_old;

% time step
dt = Control.dt;

%% Matrix partitioning
% matrices time discretization
Kpubar = Kpu./dt;
Kppbar = Kpp + S./dt;
Kupbar = -Kup;
Kpnbar = Kpn./dt;
Knubar = Knu./dt;
Knnbar = Knn./dt;

% matrix partitioning
Kuu_EE = Kuu(BC.fixed_u, BC.fixed_u);
Kuu_EF = Kuu(BC.fixed_u, BC.free_u);
Kuu_FE = Kuu(BC.free_u, BC.fixed_u);
Kuu_FF = Kuu(BC.free_u, BC.free_u);

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

Kun_EE = Kun(BC.fixed_u, BC.fixed_n);
Kun_EF = Kun(BC.fixed_u, BC.free_n);
Kun_FE = Kun(BC.free_u, BC.fixed_n);
Kun_FF = Kun(BC.free_u, BC.free_n);

Kpn_EE = Kpnbar(BC.fixed_p, BC.fixed_n);
Kpn_EF = Kpnbar(BC.fixed_p, BC.free_n);
Kpn_FE = Kpnbar(BC.free_p, BC.fixed_n);
Kpn_FF = Kpnbar(BC.free_p, BC.free_n);

Knu_EE = Knubar(BC.fixed_n, BC.fixed_u);
Knu_EF = Knubar(BC.fixed_n, BC.free_u);
Knu_FE = Knubar(BC.free_n, BC.fixed_u);
Knu_FF = Knubar(BC.free_n, BC.free_u);

Knp_EE = Knp(BC.fixed_n, BC.fixed_p);
Knp_EF = Knp(BC.fixed_n, BC.free_p);
Knp_FE = Knp(BC.free_n, BC.fixed_p);
Knp_FF = Knp(BC.free_n, BC.free_p);

Knn_EE = Knnbar(BC.fixed_n, BC.fixed_n);
Knn_EF = Knnbar(BC.fixed_n, BC.free_n);
Knn_FE = Knnbar(BC.free_n, BC.fixed_n);
Knn_FF = Knnbar(BC.free_n, BC.free_n);

% matrices reassemble
% KEE = [Kuu_EE, Kup_EE;
%     Kpu_EE, Kpp_EE];
% KEF = [Kuu_EF, Kup_EF;
%     Kpu_EF, Kpp_EF];
% KFE = [Kuu_FE, Kup_FE;
%     Kpu_FE, Kpp_FE];
% KFF = [Kuu_FF, Kup_FF;
%     Kpu_FF, Kpp_FF];

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
fpbar = fp + Kpu*u_old./dt + S*p_old./dt + Kpn*n_old./dt;

fnbar = Knu*u_old./dt + Knn*n_old./dt;

% partitioning vectors
fuF = fu(BC.free_u);
fpF = fpbar(BC.free_p);
fnF = fnbar(BC.free_n);

uE = BC.fixed_u_value;
pE = BC.fixed_p_value;
nE = BC.fixed_n_value;

% dE = [uE;pE];
% fF = [fuF; fpF];

dE = [uE;pE;nE];
fF = [fuF; fpF; fnF];

%% Solve linear system
% solve for displacement and pressure
dF = KFF\(fF - KFE *dE);

% solve for reactions
rE = KEE*dE + KEF*dF;

%% Store variables
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
%%% new
n(BC.fixed_n, 1) = nE;
n(BC.free_n, 1) = nF;
%%% end

%% Gradients
udot = (u - Iteration.u_old)/dt;
pdot = (p - Iteration.p_old)/dt;
%%% new
ndot = (n - Iteration.n_old)/dt;
end