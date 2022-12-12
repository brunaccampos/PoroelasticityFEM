function [u, udot, p, pdot, reactions] = SolverSteady(Kuu, Kup, Kpp, S, Mesh, BC, Control, Iteration)
% Solve linear system for quasi-steady case
% Input parameters: coupled matrices, BC, Control, Iteration
% Outputs:
%   u: solid displacement
%   udot: solid velocity
%   p: fluid pressure

%% Iteration data
u_old = Iteration.u_old;
p_old = Iteration.p_old;

% time step
dt = Control.dt;

%% Matrix partitioning
% matrices time discretization
Kpu = (Kup.')*dt;
Kpp = Kpp + S./dt;

% matrix partitioning
Kuu_EE = Kuu(BC.fixed_u, BC.fixed_u);
Kuu_EF = Kuu(BC.fixed_u, BC.free_u);
Kuu_FE = Kuu(BC.free_u, BC.fixed_u);
Kuu_FF = Kuu(BC.free_u, BC.free_u);

Kup_EE = Kup(BC.fixed_u, BC.fixed_p);
Kup_EF = Kup(BC.fixed_u, BC.free_p);
Kup_FE = Kup(BC.free_u, BC.fixed_p);
Kup_FF = Kup(BC.free_u, BC.free_p);

Kpu_EE = Kpu(BC.fixed_p, BC.fixed_u);
Kpu_EF = Kpu(BC.fixed_p, BC.free_u);
Kpu_FE = Kpu(BC.free_p, BC.fixed_u);
Kpu_FF = Kpu(BC.free_p, BC.free_u);

Kpp_EE = Kpp(BC.fixed_p, BC.fixed_p);
Kpp_EF = Kpp(BC.fixed_p, BC.free_p);
Kpp_FE = Kpp(BC.free_p, BC.fixed_p);
Kpp_FF = Kpp(BC.free_p, BC.free_p);

% matrices reassemble
KEE = [Kuu_EE, -Kup_EE;
    Kpu_EE, Kpp_EE];
KEF = [Kuu_EF, -Kup_EF;
    Kpu_EF, Kpp_EF];
KFE = [Kuu_FE, -Kup_FE;
    Kpu_FE, Kpp_FE];
KFF = [Kuu_FF, -Kup_FF;
    Kpu_FF, Kpp_FF];

% prescribed traction (external forces)
fu = BC.fext;
% prescribed flux
fp = zeros(Mesh.ndof_p, 1);

% auxiliar terms for external forces vector
fpbar = fp + (Kup.')*u_old/dt + S*p_old/dt;

% partitioning vectors
fuF = fu(BC.free_u);
fpF = fpbar(BC.free_p);

uE = BC.fixed_u_value;
pE = BC.fixed_p_value;

dE = [uE;pE];
fF = [fuF; fpF];

%% Solve linear system
% solve for displacement and pressure
dF = KFF\(fF - KFE *dE);

% solve for reactions
reactions = KEE*dE + KEF*dF;

%% Store variables
u = [dF(1:length(BC.free_u),1); uE];
p = [pE; dF(length(BC.free_u)+1 : length(BC.free_u) + length(BC.free_p),1)];

%% Gradients
udot = (u - u_old)/dt;
pdot = (p - p_old)/dt;
end