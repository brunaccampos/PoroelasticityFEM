function [disp_er, press_er] = PorousMediaConsolidation_check(Solution, MeshU, MeshP, Material, Control, QuadU, BC)
% Calculates the error between the FE approximation and the analytical
% solution of a porous media consolidation example in one dimension
% ------------------------------------------------------------------------
% Adapted from https://github.com/GCMLab
% ------------------------------------------------------------------------

% Calculate exact solutions
[p_an, u_an] = getAnSol_coupledComp(Control, Material, MeshU, MeshP, BC);

% Calculate the error
disp_er = norm(Solution.u - u_an);
press_er = norm(Solution.p -  p_an);

end

