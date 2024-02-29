function [Iteration, Plot] = initVariables(MeshU, MeshP, MeshN, Material, Control, BC)
% ------------------------------------------------------------------------
% Initialize variables for iteration and plot stages; store initial
% conditions
% ------------------------------------------------------------------------

%% Iteration data - time domain
% initial condition
Iteration.u_old = BC.initU; % solid displacement 
Iteration.p_old = BC.initP; % pressure
Iteration.udot_old = BC.initUdot; % solid velocity

% initialize time derivatives
Iteration.u2dot_old = zeros(MeshU.nDOF, 1); % solid acceleration
Iteration.pdot_old = zeros(MeshP.nDOF, 1); % pressure gradient

% initialize load vectors
Iteration.fu_old = zeros(MeshU.nDOF, 1); % load vector
Iteration.fp_old = zeros(MeshP.nDOF, 1); % flux vector

% initialize variables for UPN model
if contains(Control.PMmodel, 'UPN')
    % initial porosity condition
    Iteration.n_old = zeros(MeshN.nDOF, 1);
    Iteration.n_old(:) = Material.n;
    Iteration.ndot_old = zeros(MeshN.nDOF, 1);
    % load vector
    Iteration.fn_old = zeros(MeshN.nDOF, 1); % flux vector
end

% initialize variables for UPU model
if contains(Control.PMmodel, 'UPU')
    Iteration.uf_old = BC.initUf; % fluid displacement
    Iteration.ufdot_old = BC.initUfdot; % fluid velocity
    Iteration.uf2dot_old = zeros(MeshU.nDOF, 1); % fluid acceleration
    Iteration.ff_old = zeros(MeshU.nDOF, 1); % load vector
end

%% Plot arrays - time domain
if isfield(Control, 'dtmin') 
    time1 = (0:Control.dtmin:Control.tlim);
    time2 = (Control.tlim + Control.dt:Control.dt:Control.tend);
    Plot.time = [time1, time2];
else
    Plot.time = (0:Control.dt:Control.tend);
end

Plot.p_time = zeros(length(Plot.time), 1); % fluid pressure
Plot.pan_time = zeros(length(Plot.time), 1); % analytic fluid pressure (quasi-steady case)
Plot.uan_time = zeros(length(Plot.time),1); % analytic solid displacement
Plot.u_time = zeros(length(Plot.time), 1); % solid displacement
Plot.udot_time = zeros(length(Plot.time), 1); % solid velocity
Plot.u2dot_time = zeros(length(Plot.time), 1); % solid acceleration

if Control.plotansol
    aux1 = Control.p_an(0);
    aux2 = Control.u_an(0);
    Plot.pan_time(1,1) = aux1(Control.plotp);
    Plot.uan_time(1,1) = aux2(Control.plotu);
end

% store initial conditions
Plot.u_time(1,1) = BC.initU(Control.plotu,1);
Plot.p_time(1,1) = BC.initP(Control.plotp,1);
Plot.udot_time(1,1) = BC.initUdot(Control.plotu,1);

if contains(Control.PMmodel, 'UPU')
    Plot.uf_time = zeros(length(Plot.time), 1); % fluid displacement
    Plot.ufdot_time = zeros(length(Plot.time), 1); % fluid velocity
    Plot.uf2dot_time = zeros(length(Plot.time), 1); % fluid acceleration
    Plot.ufan_time = zeros(length(Plot.time),1); % analytic fluid displacement
    
    % store initial conditions
    Plot.uf_time(1,1) = BC.initUf(Control.plotu,1);
    Plot.ufdot_time(1,1) = BC.initUfdot(Control.plotu,1);
    if Control.plotansol
        aux3 = Control.uf_an(0); 
        Plot.ufan_time(1,1) = aux3(Control.plotu);
    end
end

% porosity
if contains(Control.PMmodel, 'UPN')
    Plot.n_time = zeros(length(Plot.time), 1);
    Plot.n_time(1,1) = Material.n;
end  

%% Plot arrays - synthetics
Plot.u_synthetic = zeros(length(Plot.time), length(Control.ploturow)); % solid displacement
Plot.udot_synthetic = zeros(length(Plot.time), length(Control.ploturow)); % solid velocity
Plot.p_synthetic = zeros(length(Plot.time), length(Control.plotprow)); % fluid pressure
if contains(Control.PMmodel, 'UPN')
    Plot.n_synthetic = zeros(length(Plot.time), length(Control.plotprow)); % porosity
end
if contains(Control.PMmodel, 'UPU')
    Plot.uf_synthetic = zeros(length(Plot.time), length(Control.ploturow)); % fluid displacement
    Plot.ufdot_synthetic = zeros(length(Plot.time), length(Control.ploturow)); % fluid pressure
end

%% Plot arrays - space
Plot.p_space = zeros(length(Control.plotp),1);
Plot.u_space = zeros(length(Control.plotu),1);
Plot.udot_space = zeros(length(Control.plotu),1);
Plot.u2dot_space = zeros(length(Control.plotu),1);
Plot.pan_space = zeros(length(Control.plotp),1);
Plot.uan_space = zeros(length(Control.plotu),1);

if contains(Control.PMmodel, 'UPU')
    Plot.uf_space = zeros(length(Control.plotu),1);
    Plot.ufdot_space = zeros(length(Control.plotu),1);
    Plot.uf2dot_space = zeros(length(Control.plotu),1);
end

end