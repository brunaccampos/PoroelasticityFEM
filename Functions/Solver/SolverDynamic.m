function [Solution] = SolverDynamic(Kuu, Kup, Kpp, M, Mhat, S, fu, fp, BC, Control, Iteration)
% Solve linear system for dynamic case
% Input parameters: coupled matrices, BC, Control, Iteration
% Outputs:
%   u: solid displacement
%   udot: solid velocity
%   u2dot: solid acceleration
%   p: fluid pressure
%   pdot: fluid pressure gradient

%% Iteration data
if ~isempty(Iteration)
    u_old = Iteration.u_old;
    udot_old = Iteration.udot_old;
    u2dot_old = Iteration.u2dot_old;
    p_old = Iteration.p_old;
    pdot_old = Iteration.pdot_old;
else
    u_old = zeros(length(Kuu),1);
    udot_old = zeros(length(Kuu),1);
    u2dot_old = zeros(length(Kuu),1);
    p_old = zeros(length(Kpp),1);
    pdot_old = zeros(length(Kpp),1);
end

% time step
dt = Control.dt;
% time integration parameters
beta = Control.beta;
gamma = Control.gamma;
theta = Control.theta;

%% Matrix partitioning
% matrices time discretization
Kpu = Kup.';
Kuubar = Kuu + M./(beta*dt^2);
Kpubar = Kpu *gamma/(beta*dt) - Mhat/(beta*dt^2);
Kppbar = Kpp + S/(theta*dt);
Kupbar = -Kup;

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

% matrices for unknown DOFs
MFF = M(BC.free_u, BC.free_u);
KuuFF = Kuu(BC.free_u, BC.free_u);
KupFF = Kup(BC.free_u, BC.free_p);
KpuFF = Kpu(BC.free_p, BC.free_u);
KppFF = Kpp(BC.free_p, BC.free_p);
MhatFF = Mhat(BC.free_p, BC.free_u);
SFF = S(BC.free_p, BC.free_p);

% at first step: compute solid acceleration and pressure gradient
if Control.step == 1
    u2dot_old(BC.free_u) = MFF\(-KuuFF*u_old(BC.free_u) + KupFF*p_old(BC.free_p) + fu(BC.free_u));
    pdot_old(BC.free_p) = SFF\(MhatFF*u2dot_old(BC.free_u) - KpuFF*udot_old(BC.free_u) -KppFF*p_old(BC.free_p) + fp(BC.free_p));
end

% auxiliar terms for external forces vector
fuF = fu(BC.free_u) + MFF./(beta*dt^2) * u_old(BC.free_u) + MFF/(beta*dt) * udot_old(BC.free_u) + MFF * (1/(2*beta) -1) * u2dot_old(BC.free_u);
fpF = fp(BC.free_p) + (KpuFF * gamma/(beta*dt) - MhatFF /(beta*dt^2)) * u_old(BC.free_u) + (KpuFF * (gamma/beta -1) - MhatFF/(beta*dt))* udot_old(BC.free_u)  +...
    (KpuFF * dt * (gamma/(2*beta) -1) - MhatFF *(1/(2*beta)-1)) * u2dot_old(BC.free_u) + SFF * p_old(BC.free_p) /(theta *dt) + SFF * (1/theta -1) * pdot_old(BC.free_p);

fuE = fu(BC.fixed_u);
fpE = fp(BC.fixed_p);

uE = BC.fixed_u_value;
pE = BC.fixed_p_value;

dE = [uE; pE];
fE = [fuE; fpE];
fF = [fuF; fpF];

%% Solve linear system
% solve for displacement and pressure
dF = KFF\(fF - KFE *dE);

% solve for reactions
rE = KEE*dE + KEF*dF - fE;

uF = dF(1:length(BC.free_u),1);
pF = dF(length(BC.free_u)+1 : length(BC.free_u) + length(BC.free_p),1);

% system to compute residual
A = [KEE, KEF;
    KFE, KFF];
b = [rE; fF];
x = [dE; dF];

% iteration counter
it = 1;
% convergence indicator
converged = 0;

% residuals vector
NR = zeros(Control.max_it,1);

while (it <= Control.max_it) && ~converged
    %     fprintf('It. %d ', it);

    % residuals
%     Ru = Kuubar * [uE;uF] + Kupbar * [pE;pF] - [rE(1:length(BC.fixed_u));fuF];
%     Rp = Kpubar * [uE;uF] + Kppbar * [pE;pF] - [rE(length(BC.fixed_u)+1:end);fpF];
%     JEE = [Kuu_EE, Kup_EE; Kpu_EE, Kpp_EE];
%     JEF = [Kuu_EF, Kup_EF; Kpu_EF, Kpp_EF];
%     RE = -JEE * [uE;pE] - JEF * [uF;pF];
    
    % total residual
    R = A*x - b;
    RF = R(length(BC.fixed_u) + length(BC.fixed_p) + 1 : end);

    RuF = RF(1 : length(BC.free_u),1); % residual for displacement u
    RpF = RF(length(BC.free_u) + 1 : end,1); % residual for pressure p
    
    % define value to normalize u
    if ~isempty(BC.pointLoad) % use point load if any
        u_norm = abs(BC.pointLoadValue); 
    elseif ~isempty(BC.tractionNodes) && any(BC.tractionForce,'all') % use traction magnitude if any nonzero terms
        u_norm = abs(BC.traction); % value to normalize u
    else % use unit if no traction BCs
        u_norm = 1;
    end

    % define value to normalize p
    if ~isempty(BC.pointFlux) % use point flux if any
        p_norm = abs(BC.pointFluxValue);
    elseif ~isempty(BC.fluxNodes) && any (BC.fluxValue, 'all') % use flux magnitude if any nonzero terms
        p_norm = abs(BC.flux);
    else % use unit if no flux BCs
        p_norm = 1; 
    end
    
    % norm of the residuals vectors
    Ru_norm = norm(RuF)*(1/u_norm); 
    Rp_norm = norm(RpF)*(1/p_norm); 

    % total residual norm
    NR(it) = sqrt(Ru_norm^2 + Rp_norm^2); 

    % check if converged according to defined tolerance
    if NR(it) < Control.tol
        fprintf('Converged with %d it.', it);
        converged = 1;
    else
        u_norm = mean(abs(diag(Kuu_FF)));
        p_norm = mean(abs(diag(Kpp_FF)));
        KFF(1:length(BC.free_u), :) = KFF(1:length(BC.free_u), :)*(1/u_norm);
        KFF(length(BC.free_u)+1:end, :) = KFF(length(BC.free_u)+1:end, :)*(1/p_norm);
        
        RuF = RuF*(1/u_norm);
        RpF = RpF*(1/p_norm);
        RF = [RuF; RpF];

        % solve system for increment dxF
        dxF = -KFF\RF;
        % full vector dx
        dx = [dE; dxF];
        % update variables u, p
        x = x + dx;

%         figure(1)
%         semilogy(1:iter,NR(1:iter))
%         drawnow
    end

    % update iteration counter
    it = it + 1;

    % reaching maximum of iterations
    if it == Control.max_it
        disp('Convergence failed');
        break
    end
end

%% Store u/p
% xEu = x(1:length(BC.fixed_u));
% xEp = x(length(BC.fixed_u)+1 : length(BC.fixed_u)+length(BC.fixed_p));
xFu = x(length(BC.fixed_u)+length(BC.fixed_p)+1 : length(BC.fixed_u)+length(BC.fixed_p)+length(BC.free_u));
xFp = x(length(BC.fixed_u)+length(BC.fixed_p)+length(BC.free_u)+1 : end);

% force reactions
fE = rE(1:length(BC.fixed_u),1);
%flux reactions
qE = rE(length(BC.fixed_u)+1 : end,1);

u(BC.fixed_u, 1) = uE;
u(BC.free_u, 1) = xFu;
p(BC.fixed_p, 1) = pE;
p(BC.free_p, 1) = xFp;

%% Velocity and acceleration
udot = (u - u_old)*gamma/(beta*dt) - udot_old * (gamma/beta -1) - u2dot_old * dt * (gamma/(2*beta)-1);
% udot = (u - u_old)*2*gamma/(beta*dt) + udot_old * (1 - 2*gamma/beta) + u2dot_old * dt * (1 - gamma/beta);
u2dot = (u - u_old)/(beta*dt^2) - udot_old/(beta*dt) - u2dot_old * (1/(2*beta) -1);
% u2dot = (u - u_old)*2/(beta*dt^2) - udot_old*2/(beta*dt) - u2dot_old * (1/beta -1);
pdot = (p - p_old)/(theta*dt) - (1/theta - 1) * pdot_old;

%% Store variables
Solution.u = u;
Solution.udot = udot;
Solution.u2dot = u2dot;
Solution.p = p;
Solution.pdot = pdot;
Solution.fE = fE;
Solution.qE = qE;

end