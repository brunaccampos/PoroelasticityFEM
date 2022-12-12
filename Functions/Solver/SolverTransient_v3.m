function [Solution] = SolverTransient_v3(Kuu, Kup, Kpp, S, fu, fp, BC, Control, Iteration)
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
% version 3: symmetric system
% ------------------------------------------------------------------------

%% Iteration data
if ~isempty(Iteration)
    u_old = Iteration.u_old;
    p_old = Iteration.p_old;
else
    u_old = zeros(length(Kuu),1);
    p_old = zeros(length(Kpp),1);
end
% time step
dt = Control.dt;

%% Matrix partitioning
% matrices time discretization
if Control.step == 1
    Kuubar = -Kuu;
    Kpubar = Kup.';
    Kppbar = S;
    Kupbar = Kup;
else
    Kuubar = -Kuu;
    Kpubar = Kup.';
    Kppbar = Kpp * dt + S;
    Kupbar = Kup;
end

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

% auxiliar terms for external forces vector
if Control.step == 1
    fubar = - fu;
else
    fubar = Kuubar*u_old + Kupbar*p_old;
end
fpbar = fp + Kpubar*u_old + S*p_old;

% partitioning vectors
fuF = fubar(BC.free_u);
fpF = fpbar(BC.free_p);

fuE = fubar(BC.fixed_u);
fpE = fpbar(BC.fixed_p);

uE = BC.fixed_u_value;
pE = BC.fixed_p_value;

dE = [uE;pE];
fF = [fuF; fpF];
fE = [fuE; fpE];
     
%% Solve linear system
% solve for displacement and pressure

% ------------------------------------------------------------------------
% -------- Test 1: direct method
dF = KFF\(fF - KFE *dE);
% cond(full(KFF))

% -------- Test 8: ILU decomposition + bicgstabl, default ilutp
% A = KFF;
% % cond(full(A))
% b = fF - KFE*dE;
% b = sparse(b);
% setup = struct('type','ilutp');
% [L,U] = ilu(A,setup);
% dF = bicgstabl(A,b,[],[], L, U);

% ------------------------------------------------------------------------

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
u(BC.free_u, 1) = uF;
p(BC.fixed_p, 1) = pE;
p(BC.free_p, 1) = pF;

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

end