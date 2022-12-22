function [Solution] = SolverTransient_v5(Kuu, Kup, Kpp, Kpu, S, Kpn, Knn, Knu, Knp, Kun, fu, fp, fn, BC, Control, Iteration)
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
% version 5: correcting matrices for Spanos model
% ------------------------------------------------------------------------

%% Iteration data
if ~isempty(Iteration)
    u_old = Iteration.u_old;
    p_old = Iteration.p_old;
    n_old = Iteration.n_old;
else
    u_old = zeros(length(Kuu),1);
    p_old = zeros(length(Kpp),1);
    n_old = zeros(length(Knn),1);
end
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
Kunbar = Kun;

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

Knp_EE = Knp(BC.fixed_n, BC.fixed_p);
Knp_EF = Knp(BC.fixed_n, BC.free_p);
Knp_FE = Knp(BC.free_n, BC.fixed_p);
Knp_FF = Knp(BC.free_n, BC.free_p);

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
fpbar = -fp + Kpu*u_old./dt + S*p_old./dt + Kpn*n_old./dt;
fnbar = fn + Knu*u_old./dt + Knn*n_old./dt;

% partitioning vectors
fuF = fu(BC.free_u);
fpF = fpbar(BC.free_p);
fnF = fnbar(BC.free_n);

fuE = fu(BC.fixed_u);
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

%% Gradients
udot = (u - u_old)/dt;
pdot = (p - p_old)/dt;
ndot = (n - n_old)/dt;

%% Store variables
Solution.u = u;
Solution.udot = udot;
Solution.p = p;
Solution.pdot = pdot;
Solution.n = n;
Solution.ndot = ndot;
Solution.fE = fE;
Solution.qE = qE;

end