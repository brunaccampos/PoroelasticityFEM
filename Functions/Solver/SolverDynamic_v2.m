function [Solution] = SolverDynamic_v2(Kuu, Kup, Kpp, M, Mhat, S, fu, fp, BC, Control, Iteration)
% ------------------------------------------------------------------------
% Solve linear system for dynamic case
% ------------------------------------------------------------------------
% Input parameters: coupled matrices, BC, Control, Iteration
% Outputs:
%   u: solid displacement
%   udot: solid velocity
%   u2dot: solid acceleration
%   p: fluid pressure
%   pdot: fluid pressure gradient
% ------------------------------------------------------------------------
% version 2: changing the order of the matrix partitioning. Not necessary
% ------------------------------------------------------------------------

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

% if Control.step == 1
%     u2dot_old = M\(-Kuu*u_old + Kup*p_old + fu);
% end

%% Matrix partitioning
Kpu = Kup.';

KuuEE = Kuu(BC.fixed_u, BC.fixed_u);
KuuEF = Kuu(BC.fixed_u, BC.free_u);
KuuFE = Kuu(BC.free_u, BC.fixed_u);
KuuFF = Kuu(BC.free_u, BC.free_u);

KupEE = Kup(BC.fixed_u, BC.fixed_p);
KupEF = Kup(BC.fixed_u, BC.free_p);
KupFE = Kup(BC.free_u, BC.fixed_p);
KupFF = Kup(BC.free_u, BC.free_p);

KpuEE = Kpu(BC.fixed_p, BC.fixed_u);
KpuEF = Kpu(BC.fixed_p, BC.free_u);
KpuFE = Kpu(BC.free_p, BC.fixed_u);
KpuFF = Kpu(BC.free_p, BC.free_u);

KppEE = Kpp(BC.fixed_p, BC.fixed_p);
KppEF = Kpp(BC.fixed_p, BC.free_p);
KppFE = Kpp(BC.free_p, BC.fixed_p);
KppFF = Kpp(BC.free_p, BC.free_p);

MEE = M(BC.fixed_u, BC.fixed_u);
MEF = M(BC.fixed_u, BC.free_u);
MFE = M(BC.free_u, BC.fixed_u);
MFF = M(BC.free_u, BC.free_u);

SEE = S(BC.fixed_p, BC.fixed_p);
SEF = S(BC.fixed_p, BC.free_p);
SFE = S(BC.free_p, BC.fixed_p);
SFF = S(BC.free_p, BC.free_p);

MhatEE = Mhat(BC.fixed_p, BC.fixed_u);
MhatEF = Mhat(BC.fixed_p, BC.free_u);
MhatFE = Mhat(BC.free_p, BC.fixed_u);
MhatFF = Mhat(BC.free_p, BC.free_u);


KuubarEE = KuuEE + MEE./(Control.beta*(dt^2));
KpubarEE = (KpuEE).*Control.gamma/(Control.beta*dt) - MhatEE/(Control.beta * (dt^2));
KppbarEE = KppEE + SEE/(Control.theta*dt);
KupbarEE = -KupEE;

KuubarEF = KuuEF + MEF./(Control.beta*(dt^2));
KpubarEF = (KpuEF).*Control.gamma/(Control.beta*dt) - MhatEF/(Control.beta * (dt^2));
KppbarEF = KppEF + SEF/(Control.theta*dt);
KupbarEF = -KupEF;

KuubarFE = KuuFE + MFE./(Control.beta*(dt^2));
KpubarFE = (KpuFE).*Control.gamma/(Control.beta*dt) - MhatFE/(Control.beta * (dt^2));
KppbarFE = KppFE + SFE/(Control.theta*dt);
KupbarFE = -KupFE;

KuubarFF = KuuFF + MFF./(Control.beta*(dt^2));
KpubarFF = (KpuFF).*Control.gamma/(Control.beta*dt) - MhatFF/(Control.beta * (dt^2));
KppbarFF = KppFF + SFF/(Control.theta*dt);
KupbarFF = -KupFF;

KEE = [KuubarEE, KupbarEE;
    KpubarEE, KppbarEE];
KEF = [KuubarEF, KupbarEF;
    KpubarEF, KppbarEF];
KFE = [KuubarFE, KupbarFE;
    KpubarFE, KppbarFE];
KFF = [KuubarFF, KupbarFF;
    KpubarFF, KppbarFF];


% compute acceleration at first step
if Control.step == 1
    u2dot_old(BC.free_u) = M(BC.free_u, BC.free_u)\(-KuuFF*u_old(BC.free_u) + KupFF*p_old(BC.free_p) + fu(BC.free_u));
end

% auxiliar terms for external forces vector
fubar = fu + M./(Control.beta*(dt^2)) * u_old + M/(Control.beta*dt) * udot_old + M * (1/(2*Control.beta) -1) * u2dot_old;
fpbar = fp + ((Kup.') * Control.gamma/(Control.beta*dt) - Mhat/(Control.beta*(dt^2))) * u_old + ((Kup.') * (Control.gamma/Control.beta -1) - Mhat/(Control.beta*dt))* udot_old  +...
    (Kup' * dt * (Control.gamma/(2*Control.beta) -1) - Mhat *(1/(2*Control.beta)-1)) * u2dot_old + S*p_old/(Control.theta *dt) + S*(1/Control.theta -1) * pdot_old;

% partitioning vectors
% fuF = fubar(BC.free_u);
% fpF = fpbar(BC.free_p);


fuF = fu(BC.free_u) + MFF./(Control.beta*(dt^2)) * u_old(BC.free_u) + MFF/(Control.beta*dt) * udot_old(BC.free_u) + MFF * (1/(2*Control.beta) -1) * u2dot_old(BC.free_u);
fpF = fp(BC.free_p) + (KpuFF * Control.gamma/(Control.beta*dt) - MhatFF/(Control.beta*(dt^2))) * u_old(BC.free_u) + (KpuFF * (Control.gamma/Control.beta -1) - MhatFF/(Control.beta*dt))* udot_old(BC.free_u)  +...
    (KpuFF * dt * (Control.gamma/(2*Control.beta) -1) - MhatFF *(1/(2*Control.beta)-1)) * u2dot_old(BC.free_u) + SFF * p_old(BC.free_p)/(Control.theta *dt) + SFF * (1/Control.theta -1) * pdot_old(BC.free_p);


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

%% Store u/p
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
udot = (u - u_old)*Control.gamma/(Control.beta*Control.dt) - udot_old * (Control.gamma/Control.beta -1) - u2dot_old * Control.dt * (Control.gamma/(2*Control.beta)-1);
u2dot = (u - u_old)/(Control.beta*Control.dt^2) - udot_old/(Control.beta*Control.dt) - u2dot_old * (1/(2*Control.beta) -1);
pdot = (p - p_old)/(Control.theta*Control.dt) - (1/Control.theta - 1) * pdot_old;

%% Store variables
Solution.u = u;
Solution.udot = udot;
Solution.u2dot = u2dot;
Solution.p = p;
Solution.pdot = pdot;
Solution.fE = fE;
Solution.qE = qE;

end