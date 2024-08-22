% Porous Media Simulation
% ------------------------------------------------------------------------
% Created by Bruna Campos
% bccampos@uwaterloo.ca
% Department of Civil Engineering, University of Waterloo
% January 2022
% ------------------------------------------------------------------------
% Reference: https://github.com/GCMLab
% ------------------------------------------------------------------------

%% Delete past files
if plot2vtk
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
[Material, MeshU, MeshP, MeshN, BC, Control] = feval(config_name, ConfigDir, progress_on, meshfilename, nelements);

%% Set defaults
[Material, BC, Control] = setDefaults(Material, MeshU, MeshP, BC, Control);

%% Quadrature points
QuadU = GlobalQuad(MeshU, Control);
QuadP = GlobalQuad(MeshP, Control);

%% Solution parameters
Control.t = 0;  % initial simulation time
Control.step = 0; % step counter

%% Initialize variables
vars_type = [extractAfter(Control.PMmodel, 'Spanos_'), extractAfter(Control.PMmodel, 'Biot_')]; 
initVars_function = append('initVariables', vars_type);
[Iteration, Plot] = feval(initVars_function, MeshU, MeshP, MeshN, Material, Control, BC);

%% Select model and evaluate main function
[~,modeltype] = fileparts(Control.PMmodel);
main_type = append('main',modeltype);
if Control.freqDomain
    main_type = append(main_type,'_Freq');
end
feval(main_type);

%% Compute stress - flux - porosity
% compute solid strain and stress
[strain, stress] = ComputeSolidStress(Material, MeshU, Solution.u);
Solution.strain = strain;
Solution.stress = stress;

% compute pressure gradient and flux
[gradp, flux] = ComputeFluidFlux(Material, MeshP, Solution.p);
Solution.gradp = gradp;
Solution.flux = flux;

% compute flux strain and stress
if contains(Control.PMmodel, 'UPU')
    [strainf, stressf] = ComputeSolidStress(Material, MeshU, Solution.uf);
    Solution.strainf = strainf;
    Solution.stressf = stressf;
end

% compute porosity
if contains(Control.PMmodel, 'UPU') || contains(Control.PMmodel, 'UPV') || contains(Control.PMmodel, 'UPW')
    [eta, etadot] = ComputePorosity(Material, MeshU, Solution, Control);
    Solution.n = eta;
    Solution.ndot = etadot;
end

%% Post processing
% plot figures
if plotGraphs_on && ~Control.freqDomain
    PlotGraphs(Solution, [], Material, MeshU, MeshP, MeshN, Control, Plot, saveGraphs_on);
elseif plotGraphs_on && Control.freqDomain
    PlotGraphs(Solution, SolutionFreq, Material, MeshU, MeshP, MeshN, Control, Plot, saveGraphs_on);
end

% plot synthetics
if Control.plotSyntheticsON
    PlotSynthetics(MeshU, MeshP, MeshN, Plot, Control);
end
    
% plot natural frequencies
if Control.freqDomain
    PlotModeShapes(phi_u, omega2_u, phi_p, omega2_p, MeshU, MeshP, Control, BC, config_name, vtk_dir);
end

% compute error
if saveMatData_on && Control.plotansol      
    % compute strain and stress (analytical)
    [strain_an, stress_an] = ComputeSolidStress(Material, MeshU, Plot.uan_space);
    Solution.strain_an = strain_an;
    Solution.stress_an = stress_an;

    % compute flux (analytical)
    [gradp_an, flux_an] = ComputeFluidFlux(Material, MeshP, Plot.pan_space);
    Solution.gradp_an = gradp_an;
    Solution.flux_an = flux_an;
   
    if contains(Control.PMmodel, 'UPU')
        % compute strain and stress (analytical)
        [strainf_an, stressf_an] = ComputeSolidStress(Material, MeshU, Plot.ufan_space);
        Solution.strainf_an = strainf_an;
        Solution.stressf_an = stressf_an;
    end

    % error
    if contains(Control.PMmodel, 'UPU')
        [ErrorComp] = ComputeMeshSizeError_UPU(MeshU, MeshP, Solution, Plot, Control);
    else
        [ErrorComp] = ComputeMeshSizeError_UP(MeshU, MeshP, Solution, Plot, Control);
    end

    % save results
    save('Results.mat', 'ErrorComp', 'Control', 'Plot');
end
