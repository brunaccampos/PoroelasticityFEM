% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function [p_an, u_an, uf_an] = getAnSol_uncoupled_UPU(Control, ~, ~, ~, ~)
% Return analytical solution for uncoupled problem (e.g. elasticity, heat
% transfer

p_an = Control.p_an(Control.t);
u_an = Control.u_an(Control.t);
uf_an = Control.uf_an(Control.t);

end