% ------------------------------------------------------------------------
% Runs unit Test 12 - Poroelasticity consolidation as a part of RunTests
% ------------------------------------------------------------------------
% Verify the FEA solution of a one-dimensional consolidation problem in
% porous media theory by comparing it to its analytical solution
% ------------------------------------------------------------------------
% Adapted from https://github.com/GCMLab (Acknowledgements: Bruce Gee)
% ------------------------------------------------------------------------

% Create test VTK folder
if plot2vtk
    vtk_dir = fullfile(VTKFolder,'\Test12');
    if ~isfolder(vtk_dir)
        mkdir(vtk_dir)
    end
end

fprintf('\n\n Test 12: Porous Media - consolidation\n')

%% Step 1 - Run Simulation
config_name = 'PorousMedia_Consolidation';
meshfilename = '';
nelements = 0;
main

%% Step 2 - Check results
[disp_er, press_er] = PorousMediaConsolidation_check(Solution, MeshU, MeshP, Material, Control, QuadU, BC);

fprintf('\nPM Consolidation: Displacement error is %.2f',disp_er);
fprintf('\nPM Consolidation: Pressure error is %.2f',press_er);

convergence_tolerance = 1e-10;
if disp_er <= convergence_tolerance && press_er <= convergence_tolerance
    test_pass = 1;
else
    test_pass = 0;
end

%% Step 3 - Output results
if test_pass
    fprintf('\nPASS Porous Media Consolidation\n')
else
    fprintf('\n\nFAIL Porous Media Consolidation\n')
end
testpasssummary(12) = test_pass;

%% Step 4 - Cleanup
clearvars -except  curDir  ConfigDir ...
    ntests testpasssummary...
    plot2vtk VTKFolder progress_on...
    saveGraphs_on saveMatData_on plotGraphs_on saveVideo_on