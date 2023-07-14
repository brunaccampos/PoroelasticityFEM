function [Material, BC, Control] = setDefaults(Material, BC, Control)
% Set defaults

%% Material model
if ~isfield(Material, 'constLaw')
    Material.constLaw = 'PlaneStress';
end

if ~isfield(Material, 'lumpedMass')
    Material.lumpedMass = 0;
end

if ~isfield(Material, 'lumpedDamping')
    Material.lumpedDamping = 0;
end

%% Boundary/Initial conditions
if ~isfield(BC, 'initU')
    BC.initU = [];
end

if ~isfield(BC, 'initP')
    BC.initP = [];
end

%% Control options
if ~isfield(Control, 'rampLoad')
   Control.rampLoad = 0; 
end

end