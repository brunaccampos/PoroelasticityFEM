% Porous Media Simulation
% ------------------------------------------------------------------------
% Created by Bruna Campos
% bccampos@uwaterloo.ca
% Department of Civil Engineering, University of Waterloo
% January 2022
% ------------------------------------------------------------------------
% Reference: https://github.com/GCMLab
% ------------------------------------------------------------------------
% version 5 (12/12/22): creating new main functions
% ------------------------------------------------------------------------

%% Delete past files
if plot2vtk || plot2csv_on
    if progress_on
        fprintf('%.2f: Deleting old vtk and csv files...\n', toc);
    end
    sol = fullfile(vtk_dir, '*');
    delete(sol)
end

%% Input problem
if progress_on
    fprintf('%.2f: Reading config file...\n', toc);
end
[Material, MeshU, MeshP, MeshN, BC, Control] = feval(config_name, ConfigDir, progress_on);

%% Quadrature points
QuadU = GlobalQuad(MeshU, Control);
QuadP = GlobalQuad(MeshP, Control);

%% Solution parameters
Control.t = 0;  % initial simulation time
Control.step = 0; % step counter

%% Select model name and type
if Control.Biotmodel
    if Control.steady
        main_BiotTransient;
    else
        main_BiotDynamic;
    end
else
    if Control.steady
        main_SpanosTransient;
    else
        main_SpanosDynamic;
    end
end

%% Post processing
% plot figures
if Control.freqDomain
    PlotGraphs(Solution, SolutionFreq, Material, MeshU, MeshP, MeshN, Control, Plot, saveGraphs_on);
else
    PlotGraphs(Solution, [], Material, MeshU, MeshP, MeshN, Control, Plot, saveGraphs_on);
end

% plot natural frequencies
if Control.freqDomain
    PlotModeShapes(phi_u, omega2_u, phi_p, omega2_p, MeshU, MeshP, Control, BC, config_name, vtk_dir);
end

% export CSV file
if plot2csv_on
    Plot2csv(config_name, vtk_dir, Plot);
end

% compute error
if saveMatData_on && Control.plotansol
    [ErrorComp] = ComputeMeshSizeError(Material, BC, Control, MeshU, MeshP, Solution, Plot);
    save('Results.mat', 'ErrorComp');
end

% store results for error computation
if saveMatData_on && Control.plotansol
    u = Solution.u;
    p = Solution.p;
    q = ComputeFluidFlux(Material, MeshP, p);
    u_an = Plot.uan_space;
    p_an = Plot.pan_space;
    q_an = ComputeFluidFlux(Material, MeshP, p_an);
    save('Results.mat', 'Material', 'Control', 'Plot', 'QuadU', 'QuadP', 'MeshU', 'u', 'u_an', 'MeshP', 'p', 'p_an', 'q', 'q_an');
end
