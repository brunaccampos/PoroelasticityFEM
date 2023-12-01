function [Iteration, Plot] = initVariables_Freq(phi_u, phi_p, phi_n, MeshU, MeshP, MeshN, Material, Control, BC, M)
% ------------------------------------------------------------------------
% Initialize variables for iteration and plot stages; store initial
% conditions
% ------------------------------------------------------------------------

%% Plot arrays - time domain
if isfield(Control, 'dtmin') 
    time1 = (0:Control.dtmin:Control.tlim);
    time2 = (Control.tlim + Control.dt:Control.dt:Control.tend);
    Plot.time = [time1, time2];
else
    Plot.time = (0:Control.dt:Control.tend);
end

%% Iteration data - frequency domain
Iteration.uF_old = zeros(MeshU.nDOF, 1); % displacement variable storage
Iteration.uFdot_old = zeros(MeshU.nDOF, 1); % solid velocity
Iteration.pF_old = zeros(MeshP.nDOF, 1); % pressure variable storage
Iteration.uF2dot_old = zeros(MeshU.nDOF, 1); % solid acceleration
Iteration.pFdot_old = zeros(MeshP.nDOF, 1); % pressure gradient

% initial conditions for generalized displacements
Iteration.xuF_old = (phi_u) \ Iteration.uF_old(BC.free_u);
Iteration.xuFdot_old = (phi_u) \ Iteration.uFdot_old(BC.free_u);
Iteration.xuF2dot_old = zeros(length(BC.free_u),1);
Iteration.xpF_old = (phi_p) \ Iteration.pF_old(BC.free_p);
Iteration.xpFdot_old = (phi_p) \ Iteration.pFdot_old(BC.free_p);
if contains(Control.PMmodel, 'UPN')
    Iteration.nF_old = zeros(MeshN.nDOF, 1); % porosity variable storage
    Iteration.nF_old(:) = Material.n;
    Iteration.nFdot_old = zeros(MeshN.nDOF, 1); % porosity gradient
    Iteration.xnF_old = (phi_n) \ Iteration.nF_old(BC.free_n);
    Iteration.xnFdot_old = (phi_n) \ Iteration.nFdot_old(BC.free_n);
end

if nargin > 9
    Iteration.xuF_old = (phi_u.') * M(BC.free_u, BC.free_u) * Iteration.uF_old(BC.free_u);
    Iteration.xuFdot_old = (phi_u.') * M(BC.free_u, BC.free_u) * Iteration.uFdot_old(BC.free_u);
end

%% Plot arrays - frequency domain
Plot.pF = zeros(length(Plot.time), 1); % fluid pressure
Plot.uF = zeros(length(Plot.time), 1); % solid displacement

% store initial conditions
Plot.uF(1,1) = BC.initU(Control.plotu,1);
Plot.pF(1,1) = BC.initP(Control.plotp,1);

if contains(Control.PMmodel, 'UPU')
    Plot.uFdot = zeros(length(Plot.time), 1); % solid velocity
    Plot.ufF(1,1) = BC.initUf(Control.plotu,1);
end

end