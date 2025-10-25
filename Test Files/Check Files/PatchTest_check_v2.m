% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function [press_er, flux_er] = PatchTest_check_v2(d, flux, Mesh, BC, Material)
%   Calculates the error related to pressure and flux for using a patch test
%   p = -q/k* (x+y)
%   applied flux q in both x and y directions
%   flux = [q*ones(1,Mesh.nn); q*ones(1,Mesh.nn); zeros(1,Mesh.nn)]

% Calculate exact solutions
q = BC.flux;
x = Mesh.coords(:,1);
y = Mesh.coords(:,2);
d_exact = -BC.flux/Material.M(1).kf.* (x + y);
flux_exact = [q*ones(1,Mesh.nn); q*ones(1,Mesh.nn)];

% Calculate the error
press_er = norm(d - d_exact);
flux_er = sqrt(sum((flux.' - flux_exact).^2,'all'));

end

