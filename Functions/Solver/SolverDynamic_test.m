function [Solution] = SolverDynamic_test(Kuu, Kup, Kpp, M, Mhat, S, fu, fp, BC, Control, Iteration)
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
% version 5: adapt for velocity BC in elasticity problem
% ------------------------------------------------------------------------

%% Iteration data
if ~isempty(Iteration)
    u_old = Iteration.u_old;
    udot_old = Iteration.udot_old;
    u2dot_old = Iteration.u2dot_old;
else
    u_old = zeros(length(Kuu),1);
    udot_old = zeros(length(Kuu),1);
    u2dot_old = zeros(length(Kuu),1);
end

u2dot = zeros(length(u_old),1);
udot = zeros(length(u_old),1);
u = zeros(length(u_old),1);
% time step
dt = Control.dt;

% time integration parameters
beta = Control.beta;
gamma = Control.gamma;

%% Matrix partitioning
% matrix partitioning
KEE = Kuu(BC.fixed_u, BC.fixed_u);
KEF = Kuu(BC.fixed_u, BC.free_u);
KFE = Kuu(BC.free_u, BC.fixed_u);
KFF = Kuu(BC.free_u, BC.free_u);

MEE = M(BC.fixed_u, BC.fixed_u);
MEF = M(BC.fixed_u, BC.free_u);
MFE = M(BC.free_u, BC.fixed_u);
MFF = M(BC.free_u, BC.free_u);

% matrices time discretization
KFFbar = KFF + MFF./(beta*dt^2);

% boundary condition
u_old(BC.fixed_u) = BC.fixed_u_value * Control.t;
udot_old(1,1) = 1;
udot(1,1) = 1;

u2dot(BC.fixed_u) = (udot(BC.fixed_u) - udot_old(BC.fixed_u))./(gamma*dt) + (1 - 1/gamma)*u2dot_old(BC.fixed_u);
u(BC.fixed_u) = u_old(BC.fixed_u) + (1 - beta/gamma)*dt * udot_old(BC.fixed_u) + beta*dt/gamma * udot(BC.fixed_u) + (1/2 - beta/gamma)*dt^2 * u2dot_old(BC.fixed_u);

% auxiliar terms for external forces vector
fuF = fu(BC.free_u) + MFF * (u_old(BC.free_u)./(beta*dt^2) + udot_old(BC.free_u)./(beta*dt) + u2dot_old(BC.free_u) * (1/(2*beta) -1)) - ...
    KFE * u(BC.fixed_u) - MFE * u2dot(BC.fixed_u);

fuE = fu(BC.fixed_u);

%% Solve linear system
% solve for displacement and pressure

% ------------------------------------------------------------------------
% -------- Test 1: direct method
dF = KFFbar\fuF;

% -------- Test 8: ILU decomposition + bicgstabl, default ilutp
% A = KFF;
% % cond(full(A))
% b = fF - KFE*dE;
% b = sparse(b);
% setup = struct('type','ilutp');
% [L,U] = ilu(A,setup);
% dF = bicgstabl(A,b,[],[], L, U);

% ------------------------------------------------------------------------

%% Store u/p
% u(BC.fixed_u, 1) = u;
% u(BC.free_u, 1) = dF;

%% Velocity and acceleration
u2dot(BC.free_u) = (u(BC.free_u) - u_old(BC.free_u))/(beta*dt^2) - udot_old(BC.free_u)/(beta*dt) - u2dot_old(BC.free_u) * (1/(2*beta) -1);
udot(BC.free_u) = udot_old(BC.free_u) + dt * (1-gamma) * u2dot_old(BC.free_u) + gamma * dt * u2dot(BC.free_u);

%% Store variables
Solution.u = u;
Solution.udot = udot;
Solution.u2dot = u2dot;

p = zeros(length(Kpp),1);
Solution.p = p;
pdot = zeros(length(Kpp),1);
Solution.pdot = pdot;
end