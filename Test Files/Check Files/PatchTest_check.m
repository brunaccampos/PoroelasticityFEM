% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function [disp_er, stress_er, reaction_er] = PatchTest_check(d, stress, Fext, Mesh, BC, Material)
%   Calculates the error related to displacements, reaction forces, and
%   stress for using a patch test
%   ux = (1-nu)*t/E*x
%   uy = (1-nu)*t/E*y
%   applied traction of t in both x and y directions
%   stress = [t*ones(1,Mesh.nn); t*ones(1,Mesh.nn); zeros(1,Mesh.nn)]

% Calculate exact solutions
sigma = BC.traction;
x = Mesh.coords(:,1);
y = Mesh.coords(:,2);
d_exact = zeros(2*Mesh.nn,1);
d_exact(1:2:end) = (1-Material.M(1).nu)*BC.traction/Material.M(1).E.*x;
d_exact(2:2:end) = (1-Material.M(1).nu)*BC.traction/Material.M(1).E.*y;
stress_exact = [sigma*ones(1,Mesh.nn); sigma*ones(1,Mesh.nn); zeros(1,Mesh.nn)];

% Calculate the error
disp_er = norm(d - d_exact);
stress_er = sqrt(sum((stress.' - stress_exact).^2,'all'));
reaction_er = abs(sum(Fext));

end

