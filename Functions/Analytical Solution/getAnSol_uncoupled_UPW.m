function [p_an, u_an, w_an] = getAnSol_uncoupled_UPW(Control, ~, ~, ~, ~)
% ------------------------------------------------------------------------
% Return analytical solution for uncoupled problem (e.g. elasticity, heat
% transfer
% ------------------------------------------------------------------------

p_an = Control.p_an(Control.t);
u_an = Control.u_an(Control.t);
w_an = Control.w_an(Control.t);

end