% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function PlotGraphs(Solution, SolutionFreq, Material, MeshU, MeshP, MeshN, Control, Plot, saveGraphs_on)
% Plot graphic results

if contains(Control.PMmodel, 'Tr')
    switch MeshU.nsd
        case 1
            % One-dimensional transient
            PlotGraphs1D_Tr(Solution, SolutionFreq, Material, MeshU, MeshP, MeshN, Control, Plot, saveGraphs_on);
        case 2
            % Two-dimensional transient
            PlotGraphs2D_Tr(MeshU, MeshP, MeshN, Control, Plot, Material, saveGraphs_on);
    end
else
    switch MeshU.nsd
        case 1
            % One-dimensional dynamic
            PlotGraphs1D_Dyn(Solution, SolutionFreq, MeshU, MeshP, MeshN, Control, Plot, Material, saveGraphs_on);
        case 2
            % Two-dimensional dynamic
            PlotGraphs2D_Dyn(MeshU, MeshP, MeshN, Control, Plot, Material, saveGraphs_on);
    end
end

end
