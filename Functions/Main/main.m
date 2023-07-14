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

%% Set defaults
[Material, BC, Control] = setDefaults(Material, BC, Control);

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

% plot synthetics
if isfield(Control, 'depthplot') && contains(Control.PMmodel, 'Dyn')
    PlotSynthetics(MeshU, MeshP, MeshN, Plot, Control);
end
    
% plot natural frequencies
if Control.freqDomain
    PlotModeShapes(phi_u, omega2_u, phi_p, omega2_p, MeshU, MeshP, Control, BC, config_name, vtk_dir);
end

% compute error
if saveMatData_on && Control.plotansol
    % symbolic analytical results
    if ~Control.uncoupled
        if any(Material.Minv)
            [pan_symb, uan_symb] = getAnalyticResult_Comp_Symb(Material, MeshU, BC, Control);
        else
            [pan_symb, uan_symb] = getAnalyticResult_Incomp_Symb(Material, MeshU, BC, Control);
        end
        Control.uan_symb = uan_symb;
        Control.pan_symb = pan_symb;
    end

    % compute strain and stress
    [strain, stress] = ComputeSolidStress(Material, MeshU, Solution.u);
    [strain_an, stress_an] = ComputeSolidStress(Material, MeshU, Plot.uan_space);
    % store results
    Solution.strain = strain;
    Solution.stress = stress;
    Solution.strain_an = strain_an;
    Solution.stress_an = stress_an;

    % compute flux
    [gradp, flux] = ComputeFluidFlux(Material, MeshP, Solution.p);
    [gradp_an, flux_an] = ComputeFluidFlux(Material, MeshP, Plot.pan_space);
    % store results
    Solution.gradp = gradp;
    Solution.flux = flux;
    Solution.gradp_an = gradp_an;
    Solution.flux_an = flux_an;
   
    % error
    [ErrorComp] = ComputeMeshSizeError(MeshU, MeshP, Solution, Plot, Control);
    % save results
    save('Results.mat', 'ErrorComp');
end
