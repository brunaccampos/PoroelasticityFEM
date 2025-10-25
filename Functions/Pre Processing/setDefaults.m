% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function [Material, MeshU, MeshP, MeshN, BC, Control] = setDefaults(Material, MeshU, MeshP, MeshN, BC, Control)
% Set defaults for Material, BC, and Control structures

%% Output information
err_message = sprintf('-------------------------------------------------\n');
err_message = sprintf('%s\tfrom\t: setDefaults.m \n', err_message);
war_message = sprintf('-------------------------------------------------\n');
war_message = sprintf('%s\tfrom\t: setDefaults.m \n', war_message);

err_count = 0; % counts number of errors (missing data)
war_count = 0; % counts number of warnings (shows when a default parameter has been set)

%% Material model
err_mat = sprintf('\t\tMaterial \n');
war_mat = sprintf('\t\tMaterial \n');

% constitutive law: default plane stress
if ~isfield(Material, 'constLaw')
    war_count = war_count+1;
    war_mat = sprintf('%s\t\t\tWarning #%d\t:\t Constitutive law not defined. Set to plane stress \n',war_mat,war_count);
    Material.constLaw = 'PlaneStress';
end

% lumped mass matrix: default false
if ~isfield(Material, 'lumpedMass')
    Material.lumpedMass = 0;
end

% lumped damping matrix: default false
if ~isfield(Material, 'lumpedDamping')
    Material.lumpedDamping = 0;
end

% added mass coefficient: default zero
if ~isfield(Material.M, 'rho12')
    [Material.M(:).rho12] = deal(0);
end

% BT high-frequency correction factor: default F_BT = 1
if ~isfield(Material.M, 'F_BT')
    [Material.M(:).F_BT] = deal(1);
end

if ~isfield(MeshU, 'MatList')
    MeshU.MatList = ones(MeshU.ne, 1, 'int8');
end

if ~isfield(MeshP, 'MatList')
    MeshP.MatList = ones(MeshP.ne, 1, 'int8');
end

if contains(Control.PMmodel, 'UPN') && ~isfield(MeshN, 'MatList')
    MeshN.MatList = ones(MeshN.ne, 1, 'int8');
end

if ~isfield(MeshU, 'MatNodes')
    MeshU.MatNodes = ones(MeshU.nn,1);
end

if length(unique(MeshU.MatList))>1 && contains(Control.PMmodel, 'UPN')
    err_count = err_count+1;
    err_mat = sprintf('%s\t\t\tError #%d\t:\t This formulation does not support more than one material assigned on the domain. \n',err_mat,err_count);
end

% mapping vector
if MeshU.nsd == 1
    Material.m = 1;
elseif MeshU.nsd == 2
    Material.m = [1; 1; 0];
end

%% Boundary/Initial conditions
err_BC = sprintf('\t\tBoundary conditions \n');
war_BC = sprintf('\t\tBoundary conditions \n');

% initial condition for solid displacement field
if ~isfield(BC, 'initU')
    BC.initU = zeros(MeshU.nDOF,1);
end

% initial condition fluid pressure field
if ~isfield(BC, 'initP')
    BC.initP = zeros(MeshP.nDOF,1);
end

% initial condition for solid velocity field
if ~isfield(BC, 'initUdot')
    BC.initUdot = zeros(MeshU.nDOF,1);
end

% initial condition for fluid displacement field
if (contains(Control.PMmodel, 'UPU') || contains(Control.PMmodel, 'UPV')) && ~isfield(BC, 'initUf')
    BC.initUf = zeros(MeshU.nDOF,1);
end

% initial condition for fluid velocity field
if (contains(Control.PMmodel, 'UPU') || contains(Control.PMmodel, 'UPV')) && ~isfield(BC, 'initUfdot')
    BC.initUfdot = zeros(MeshU.nDOF,1);
end

% initial condition for relative fluid velocity field
if contains(Control.PMmodel, 'UPW') && ~isfield(BC, 'initW')
    BC.initW = zeros(MeshU.nDOF,1);
end

% boundary condition for fluid displacement
if ~isfield(BC, 'fixed_uf')
    BC.fixed_uf = [];
    BC.fixed_uf_value = @(t) zeros(length(BC.fixed_uf),1);
    BC.free_uf = setdiff(MeshU.DOF, BC.fixed_uf);
end

% traction interpolation (needed for traction applied in circular
% geometries)
if ~isfield(BC,'tractionInterp')
    BC.tractionInterp = 0;
end

% body force in fluid motion equation (upU, upv, upw model)
if ~isfield(BC,'bf') && (contains(Control.PMmodel, 'UPU') || contains(Control.PMmodel, 'UPV') || contains(Control.PMmodel, 'UPW'))
    BC.bf = @(x,t) [];
end

% body force in solid motion equation (upU, upv, upw model)
if ~isfield(BC,'bs') && (contains(Control.PMmodel, 'UPU') || contains(Control.PMmodel, 'UPV') || contains(Control.PMmodel, 'UPW'))
    BC.bs = @(x,t) [];
end

% point load  - initialize default
if ~isfield(BC, 'pointLoad')
    BC.pointLoad = @(t) [];
end

% point load - check time dependent vector
if ~isa(BC.pointLoad, 'function_handle')
    BC.pointLoad = @(t) BC.pointLoad;
end

% point flux - initialize default
if ~isfield(BC, 'pointFlux')
    BC.pointFlux = @(t) [];
end

% point flud - check time dependent vector
if ~isa(BC.pointFlux, 'function_handle')
    BC.pointFlux = @(t) BC.pointFlux;
end

if isfield(BC, 'tractionForce') && ~isa(BC.tractionForce, 'function_handle')
    BC.tractionForce = @(t) BC.tractionForce;
end

%% Control options
err_control = sprintf('\t\tControl Parameters \n');
war_control = sprintf('\t\tControl Parameters \n');

% frequency domain (mode superposition) tag: default false
if ~isfield(Control, 'freqDomain')
    Control.freqDomain = 0;
end

% ramp applied load in the beginning of simulation: default false
if ~isfield(Control, 'rampLoad')
   Control.rampLoad = 0; 
end

% plot in a row at fixed coordinate: default false
if ~isfield(Control, 'fixedDepthPlotON')
    Control.fixedDepthPlotON = 0;
end

% plot synthetics: default false
if ~isfield(Control, 'plotSyntheticsON')
    Control.plotSyntheticsON = 0;
end

% compute analytical solution: default false
if ~isfield(Control, 'plotansol')
    Control.plotansol = 0;
end

% analytical solution for fluid displacement field
if contains(Control.PMmodel, 'Tr') && ~contains(Control.PMmodel, 'UPU')
   Control.uf_an = @(t) []; 
end

% check if HHT method is used: default false
if ~isfield(Control, 'alpha')
    Control.alpha = 0;
end

% parallel processing pool
if ~isfield(Control, 'parallel')
    Control.parallel = 1;
end

% row of nodes to plot: default plot bottom nodes (DOF y) for u
if ~isfield(Control, 'ploturow') && MeshU.nsd == 1
    Control.ploturow = MeshU.DOF;
elseif ~isfield(Control, 'ploturow') && MeshU.nsd == 2
    war_count = war_count+1;
    war_control = sprintf('%s\t\t\tWarning #%d\t:\t Row of nodes to plot (u) not defined. Set as bottom y DOF \n',war_control,war_count);
    Control.ploturow = MeshU.bottom_dofy;
end

% row of nodes to plot: default plot bottom nodes for p
if ~isfield(Control, 'plotprow') && MeshU.nsd == 1
    Control.plotprow = MeshP.DOF;
elseif ~isfield(Control, 'plotprow') && MeshU.nsd == 2
    war_count = war_count+1;
    war_control = sprintf('%s\t\t\tWarning #%d\t:\t Row of nodes to plot (p) not defined. Set as bottom nodes \n',war_control,war_count);
    Control.plotprow = MeshP.bottom_nodes;
end

% node to plot: default plot node 1 in u
if ~isfield(Control,'plotu')
    Control.plotu = 1;
end

% node to plot: default plot node 1 in p
if ~isfield(Control,'plotp')
    Control.plotp = 1;
end

% scaling factor for synthetic plots
if ~isfield(Control, 'plotSyntheticScale')
    Control.plotSyntheticScale = 1;
end

%% Output error/warning
war_message = sprintf('%s %s %s %s %s', war_message, war_mat, war_BC, war_control);
err_message = sprintf('%s %s %s %s %s', err_message, err_mat, err_BC, err_control);

% check if there are any warnings
if war_count > 0
    warning('\n%s', war_message);
end

% check if there are any errors
if err_count > 0 
    error('\n%s', err_message);
end

end