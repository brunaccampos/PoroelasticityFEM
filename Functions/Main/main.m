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

%% Select model and evaluate main function
[~,modeltype] = fileparts(Control.PMmodel);
main_type = append('main',modeltype);
feval(main_type);

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
    % symbolic analytical results
    if ~Control.uncoupled
        [pan_symb, uan_symb] = getAnalyticResult_Symb(Material, MeshU, BC, Control);
        Control.uan_symb = uan_symb;
        Control.pan_symb = pan_symb;
    end

    % compute strain and stress
    [e,s] = ComputeSolidStress(Material, MeshU, Solution.u);
    [e_an,s_an] = ComputeSolidStress(Material, MeshU, Plot.uan_space);
    % store results
    Solution.e = e;
    Solution.s = s;
    Solution.e_an = e_an;
    Solution.s_an = s_an;

    % compute flux
    q = ComputeFluidFlux(Material, MeshP, Solution.p);
    q_an = ComputeFluidFlux(Material, MeshP, Plot.pan_space);
    % store results
    Solution.q = q;
    Solution.q_an = q_an;
   
    % error
    [ErrorComp] = ComputeMeshSizeError(MeshU, MeshP, Solution, Plot, Control);
    % save results
    save('Results.mat', 'ErrorComp');
end
