function [Material, BC, Control] = setDefaults(Material, BC, Control)
% ------------------------------------------------------------------------
% Set defaults for Material, BC, and Control structures
% ------------------------------------------------------------------------

%% Material model
% constitutive law: default plane stress
if ~isfield(Material, 'constLaw')
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
if ~isfield(Material, 'rho12')
    Material.rho12 = 0;
end

%% Boundary/Initial conditions
% initial condition for solid displacement field
if ~isfield(BC, 'initU')
    BC.initU = [];
end

% initial condition fluid pressure field
if ~isfield(BC, 'initP')
    BC.initP = [];
end

%% Control options
% ramp applied load in the beginning of simulation: default false
if ~isfield(Control, 'rampLoad')
   Control.rampLoad = 0; 
end

% plot in a row at fixed coordinate: default false
if ~isfield(Control, 'fixedDepthPlotON')
    Control.fixedDepthPlotON = 0;
end

end