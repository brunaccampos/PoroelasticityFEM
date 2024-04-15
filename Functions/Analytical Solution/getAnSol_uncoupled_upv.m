function [p_an, u_an, ufdot_an] = getAnSol_uncoupled_upv(Control, ~, ~, ~, ~)
% ------------------------------------------------------------------------
% Return analytical solution for uncoupled problem (e.g. elasticity, heat
% transfer
% ------------------------------------------------------------------------

p_an = Control.p_an(Control.t);
u_an = Control.u_an(Control.t);
ufdot_an = Control.ufdot_an(Control.t);

end