% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function [disp_er, press_er] = PorousMediaConsolidation_check(Solution, MeshU, MeshP, Material, Control, ~, BC)
% Calculates the error between the FE approximation and the analytical
% solution of a porous media consolidation example in one dimension

% Calculate exact solutions
[p_an, u_an] = getAnSol_coupledComp(Control, Material, MeshU, MeshP, BC);

% Calculate the error
disp_er = norm(Solution.u - u_an);
press_er = norm(Solution.p -  p_an);

end

