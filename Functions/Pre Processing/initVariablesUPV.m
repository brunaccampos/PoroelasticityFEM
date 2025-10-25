% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function [Iteration, Plot] = initVariablesUPV(MeshU, MeshP, ~, ~, Control, BC)
% Initialize variables for Iteration and Plot structures
% Model u-p-v (solid displacement, fluid pressure, fluid velocity

%% Iteration data - time domain
% initial condition
Iteration.u_old = BC.initU; % solid displacement
Iteration.p_old = BC.initP; % pressure
Iteration.udot_old = BC.initUdot; % solid velocity
Iteration.ufdot_old = BC.initUfdot; % fluid velocity

% time derivatives
Iteration.u2dot_old = zeros(MeshU.nDOF, 1); % solid acceleration
Iteration.pdot_old = zeros(MeshP.nDOF, 1); % pressure gradient
Iteration.uf2dot_old = zeros(MeshU.nDOF, 1); % fluid acceleration

% load vectors
Iteration.fu_old = zeros(MeshU.nDOF, 1); % load vector
Iteration.fp_old = zeros(MeshP.nDOF, 1); % flux vector
Iteration.ff_old = zeros(MeshU.nDOF, 1); % load vector

%% Plot arrays - time domain
if isfield(Control, 'dtmin')
    time1 = (0:Control.dtmin:Control.tlim);
    time2 = (Control.tlim + Control.dt:Control.dt:Control.tend);
    Plot.time = [time1, time2];
else
    Plot.time = (0:Control.dt:Control.tend);
end

% FE solution
Plot.p_time = zeros(length(Plot.time), 1); % fluid pressure
Plot.u_time = zeros(length(Plot.time), 1); % solid displacement
Plot.udot_time = zeros(length(Plot.time), 1); % solid velocity
Plot.u2dot_time = zeros(length(Plot.time), 1); % solid acceleration
Plot.ufdot_time = zeros(length(Plot.time), 1); % fluid velocity
Plot.uf2dot_time = zeros(length(Plot.time), 1); % fluid acceleration

% analytical solution
Plot.pan_time = zeros(length(Plot.time), 1); % analytic fluid pressure (quasi-steady case)
Plot.uan_time = zeros(length(Plot.time),1); % analytic solid displacement
Plot.ufdotan_time = zeros(length(Plot.time),1); % analytic fluid velocity

% initial conditions for FE solution
Plot.u_time(1,1) = BC.initU(Control.plotu,1);
Plot.p_time(1,1) = BC.initP(Control.plotp,1);
Plot.udot_time(1,1) = BC.initUdot(Control.plotu,1);
Plot.ufdot_time(1,1) = BC.initUfdot(Control.plotu,1);

% initial conditions for analytical solution
if Control.plotansol
    aux1 = Control.p_an(0);
    aux2 = Control.u_an(0);
    aux3 = Control.ufdot_an(0);
    Plot.pan_time(1,1) = aux1(Control.plotp);
    Plot.uan_time(1,1) = aux2(Control.plotu);
    Plot.ufdotan_time(1,1) = aux3(Control.plotu);
end

%% Plot arrays - synthetics
Plot.u_synthetic = zeros(length(Plot.time), length(Control.ploturow)); % solid displacement
Plot.udot_synthetic = zeros(length(Plot.time), length(Control.ploturow)); % solid velocity
Plot.p_synthetic = zeros(length(Plot.time), length(Control.plotprow)); % fluid pressure
Plot.ufdot_synthetic = zeros(length(Plot.time), length(Control.ploturow)); % fluid pressure

%% Plot arrays - space
Plot.p_space = zeros(length(Control.plotp),1);
Plot.u_space = zeros(length(Control.plotu),1);
Plot.udot_space = zeros(length(Control.plotu),1);
Plot.u2dot_space = zeros(length(Control.plotu),1);
Plot.pan_space = zeros(length(Control.plotp),1);
Plot.uan_space = zeros(length(Control.plotu),1);
Plot.ufdot_space = zeros(length(Control.plotu),1);
Plot.uf2dot_space = zeros(length(Control.plotu),1);
Plot.ufan_space = zeros(length(Control.plotu),1);
Plot.ufdotan_space = zeros(length(Control.plotu),1);

end