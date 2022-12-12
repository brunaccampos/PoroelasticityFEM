function [press_er, flux_er, reaction_er] = PatchTest_check_v2(d, flux, Fext, Mesh, BC, Material)
%PATCHTEST_CHECK Calculates the error
%   [disp_er, stress_er, reaction_er] = PatchTest_check(d, stress, Fext)
%   calculates the error related to displacements, reaction forces, and
%   stress for using a patch test
%
%   ----------------------------------------------------------
%   Input
%   ----------------------------------------------------------
%   d:                  Pressure vector
%   flux:               Nodal flux data
%   Fext:               External forces as the reactions
%   Mesh:               Mesh data structure
%
%   ----------------------------------------------------------
%   Output
%   ----------------------------------------------------------
%   press_er:           Error related to pressure vector
%   flux_er:            Error related to nodal flux
%   reaction_er:        Error related to reaction forces

% Patch Test
% p = -q/k* (x+y)
%   applied flux q in both x and y directions
%   flux = [q*ones(1,Mesh.nn); q*ones(1,Mesh.nn); zeros(1,Mesh.nn)]

%   ----------------------------------------------------------
% Adapted from https://github.com/GCMLab
%   ----------------------------------------------------------

% Calculate exact solutions
q = BC.flux;
x = Mesh.coords(:,1);
y = Mesh.coords(:,2);
d_exact = -BC.flux/Material.kf.* (x + y);
flux_exact = [q*ones(1,Mesh.nn); q*ones(1,Mesh.nn)];

% Calculate the error
press_er = norm(d - d_exact);
flux_er = sqrt(sum((flux.' - flux_exact).^2,'all'));
reaction_er = abs(sum(Fext));

end

