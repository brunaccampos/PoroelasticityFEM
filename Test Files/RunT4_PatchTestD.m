% ------------------------------------------------------------------------
% Runs unit Test 4 - Patch Test D as part of RunTests
% ------------------------------------------------------------------------
% For Patch Test D, all nodes are restrained and nodal pressure values
% are specfied according to the exact solution. The error between the FEA
% and exact solutions is then calculated. The FEA approximate solution
% should be exact.
% ------------------------------------------------------------------------
% Adapted from https://github.com/GCMLab (Acknowledgements: Bruce Gee)
% ------------------------------------------------------------------------

% Create test VTK folder
if plot2vtk
    vtk_dir = fullfile(VTKFolder,'\Test4');
    if ~isfolder(vtk_dir)
        mkdir(vtk_dir)
    end
end

fprintf('\n\n Test 4: Patch Test D - Q4 elements\n')

%% Step 1 - Run Simulation
config_name = 'PatchTestD';
meshfilename = 'Mesh Files\PatchTest.msh';
main

%% Step 2 - Check results
[press_er, flux_er] = PatchTest_check_v2(Solution.p, flux, MeshP, BC, Material);

fprintf('\nQ4-patch test D: Pressure error is %.2f',press_er)
fprintf('\nQ4-patch test D: Flux error is %.2f',flux_er)

convergence_tolerance = 1e-10;
if press_er <= convergence_tolerance && flux_er <= convergence_tolerance
    test_pass = 1;
else
    test_pass = 0;
end

%% Step 3 - Output results
if test_pass
    fprintf('\nPASS Patch Test D\n')
else
    fprintf('\n\nFAIL Patch Test D\n')
    return
end
testpasssummary(4) = test_pass;

%% Step 4 - Cleanup
clearvars -except  curDir  ConfigDir ...
    ntests testpasssummary...
    plot2vtk VTKFolder progress_on...
    saveGraphs_on saveMatData_on