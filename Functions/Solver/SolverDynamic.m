function [u, udot, u2dot, p, pdot, fE, qE] = SolverDynamic(Kuu, Kup, Kpp, M, Mhat, S, fu, fp, BC, Control, Iteration)
% Solve linear system for dynamic case
% Input parameters: coupled matrices, BC, Control, Iteration
% Outputs:
%   u: solid displacement
%   udot: solid velocity
%   u2dot: solid acceleration
%   p: fluid pressure
%   pdot: fluid pressure gradient

%% Iteration data
u_old = Iteration.u_old;
udot_old = Iteration.udot_old;
u2dot_old = Iteration.u2dot_old;

p_old = Iteration.p_old;
pdot_old = Iteration.pdot_old;

% time step
dt = Control.dt;

%% Matrix partitioning
% matrices time discretization
Kuu = Kuu + M/(Control.beta*dt^2);
Kpu = (Kup.')*Control.gamma/(Control.beta*dt) - Mhat/(Control.beta * dt^2);
Kpp = Kpp + S/(Control.theta*dt);
Kup = -Kup;

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
KEE = [Kuu_EE, Kup_EE;
    Kpu_EE, Kpp_EE];
KEF = [Kuu_EF, Kup_EF;
    Kpu_EF, Kpp_EF];
KFE = [Kuu_FE, Kup_FE;
    Kpu_FE, Kpp_FE];
KFF = [Kuu_FF, Kup_FF;
    Kpu_FF, Kpp_FF];

% auxiliar terms for external forces vector
fubar = fu + M/(Control.beta*dt^2) * u_old + M/(Control.beta*dt) * udot_old + M * (1/(2*Control.beta) -1) * u2dot_old;

fpbar = fp + ((-Kup.') * Control.gamma/(Control.beta*dt) - Mhat/(Control.beta*dt^2)) * u_old + (-Kup' * (Control.gamma/Control.beta -1) - Mhat/(Control.beta*dt))* udot_old  +...
    (-Kup' * dt * (Control.gamma/(2*Control.beta) -1) - Mhat *(1/(2*Control.beta)-1)) * u2dot_old + S*p_old/(Control.theta *dt) + S*(1/Control.theta -1) * pdot_old;

% partitioning vectors
fuF = fubar(BC.free_u);
fpF = fpbar(BC.free_p);

uE = BC.fixed_u_value;
pE = BC.fixed_p_value;

dE = [uE;pE];
fF = [fuF; fpF];

%% Solve linear system
% solve for displacement and pressure
dF = KFF\(fF - KFE *dE);

% solve for reactions
rE = KEE*dE + KEF*dF;

uF = dF(1:length(BC.free_u),1);
pF = dF(length(BC.free_u)+1 : length(BC.free_u) + length(BC.free_p),1);

% system to compute residual
A = [KEE, KEF;
    KFE, KFF];
b = [rE; fF];
x = [uE; pE; uF; pF];

% iteration counter
it = 1;
% convergence indicator
converged = 0;

% residuals vector
NR = zeros(Control.max_it,1);

while (it <= Control.max_it) && ~converged
    %     fprintf('It. %d ', it);

    % total residual
    R = A*x - b;
    RF = R(length(BC.fixed_u) + length(BC.fixed_p) + 1 : end);

    RuF = RF(1 : length(BC.free_u),1); % residual for displacement u
    RpF = RF(length(BC.free_u) + 1 : end,1); % residual for pressure p

    if ~isempty(BC.pointLoad)
        u_norm = abs(BC.pointLoadValue); % value to normalize u
    elseif ~isempty(BC.tractionNodes)
        u_norm = abs(BC.traction); % value to normalize u
    end

    p_norm = 1; % value to normalize p (would be zero, so use 1)

    Ru_norm = norm(RuF)*(1/u_norm); % norm of the residuals vector (u)
    Rp_norm = norm(RpF)*(1/p_norm); % norm of the residuals vector (p)

    NR(it) = sqrt(Ru_norm^2 + Rp_norm^2); % total residual norm

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

%% Store variable
xEu = x(1:length(BC.fixed_u));
xEp = x(length(BC.fixed_u)+1 : length(BC.fixed_u)+length(BC.fixed_p));
xFu = x(length(BC.fixed_u)+length(BC.fixed_p)+1 : length(BC.fixed_u)+length(BC.fixed_p)+length(BC.free_u));
xFp = x(length(BC.fixed_u)+length(BC.fixed_p)+length(BC.free_u)+1 : end);

% force reactions
fE = rE(1:length(BC.fixed_u),1);
%flux reactions
qE = rE(length(BC.fixed_u)+1 : end,1);

u(BC.fixed_u, 1) = xEu;
u(BC.free_u, 1) = xFu;
p(BC.fixed_p, 1) = xEp;
p(BC.free_p, 1) = xFp;

%% Velocity and acceleration
udot = (u - Iteration.u_old)*Control.gamma/(Control.beta*Control.dt) - Iteration.udot_old * (Control.gamma/Control.beta -1) - Iteration.u2dot_old * Control.dt * (Control.gamma/(2*Control.beta)-1);
u2dot = (u - Iteration.u_old)/(Control.beta*Control.dt^2) - Iteration.udot_old/(Control.beta*Control.dt) - Iteration.u2dot_old * (1/(2*Control.beta) -1);
pdot = (p - Iteration.p_old)/(Control.theta*Control.dt) - (1/Control.theta - 1) * Iteration.pdot_old;

end