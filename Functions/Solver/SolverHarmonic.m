function [Solution] = SolverHarmonic(Kuu, Kup, Kpp, S, M, fu, fp, BC, Control, Iteration)
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

%% Iteration data
if ~isempty(Iteration)
    u_old = Iteration.u_old;
    p_old = Iteration.p_old;
else
    % adapted for patch tests
    u_old = zeros(length(Kuu),1);
    p_old = zeros(length(Kpp),1);
end

% time step
dt = Control.dt;

%% Matrix partitioning
Kuubar = Kuu - Control.omega^2 * M;
Kupbar = -real(Kup);
Kpubar = Control.omega^2 * real(Kup.');
Kppbar = real(Kpp) - Control.omega^2 * S;

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

% matrices reassemble
KEE = [Kuu_EE, Kup_EE;
    Kpu_EE, Kpp_EE];
KEF = [Kuu_EF, Kup_EF;
    Kpu_EF, Kpp_EF];
KFE = [Kuu_FE, Kup_FE;
    Kpu_FE, Kpp_FE];
KFF = [Kuu_FF, Kup_FF;
    Kpu_FF, Kpp_FF];

% partitioning vectors
fuF = fu(BC.free_u);
fpF = fp(BC.free_p);

fuE = fu(BC.fixed_u);
fpE = fp(BC.fixed_p);

uE = BC.fixed_u_value;
pE = BC.fixed_p_value;

dE = [uE;pE];
fF = [fuF; fpF];
fE = [fuE; fpE];
     
%% Solve linear system
% solve for displacement and pressure
dF = KFF\(fF - KFE *dE);
% solve for reactions
rE = KEE*dE + KEF*dF - fE;

%% Store u/p
% unknown displacement
uF = dF(1:length(BC.free_u),1);
% unknown pressure
pF = dF(length(BC.free_u)+1 : length(BC.free_u) + length(BC.free_p),1);

% force reactions
fE = rE(1:length(BC.fixed_u),1);
%flux reactions
qE = rE(length(BC.fixed_u)+1 : end,1);

u(BC.fixed_u, 1) = uE;
u(BC.free_u, 1) = real(uF);
p(BC.fixed_p, 1) = pE;
p(BC.free_p, 1) = real(pF);

%% Gradients
udot = (u - u_old)./dt;
pdot = (p - p_old)./dt;

%% Store variables
Solution.u = u;
Solution.udot = udot;
Solution.p = p;
Solution.pdot = pdot;
Solution.fE = fE;
Solution.qE = qE;
Solution.u2dot = zeros(length(Kuu),1);
end