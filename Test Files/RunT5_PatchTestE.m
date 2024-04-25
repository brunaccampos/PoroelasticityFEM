% ------------------------------------------------------------------------
% Runs unit Test 5 - Patch Test D as part of RunTests
% ------------------------------------------------------------------------
% For Patch Test E, only nodes 1-8 (nodes in the boundaries) are restrained
% with their pressures specified according to the exact solution. The error between the FEA
% and exact solutions is then calculated. The FEA approximate solution
% should be exact.
% ------------------------------------------------------------------------
% Adapted from https://github.com/GCMLab (Acknowledgements: Bruce Gee)
% ------------------------------------------------------------------------

% Create test VTK folder
if plot2vtk
    vtk_dir = fullfile(VTKFolder,'\Test5');
    if ~isfolder(vtk_dir)
        mkdir(vtk_dir)
    end
end

fprintf('\n\n Test 5: Patch Test E - Q4 elements\n')

%% Step 1 - Run Simulation
config_name = 'PatchTestE';
meshfilename = 'Mesh Files\PatchTest.msh';
nelements = 0;
main

%% Step 2 - Check results
[press_er, flux_er] = PatchTest_check_v2(Solution.p, flux, MeshP, BC, Material);

fprintf('\nQ4-patch test E: Pressure error is %.2f',press_er)
fprintf('\nQ4-patch test E: Flux error is %.2f',flux_er)

convergence_tolerance = 1e-10;
if press_er <= convergence_tolerance && flux_er <= convergence_tolerance
    test_pass = 1;
else
    test_pass = 0;
end

%% Step 3 - Output results
if test_pass
    fprintf('\nPASS Patch Test E\n')
else
    fprintf('\n\nFAIL Patch Test E\n')
    return
end
testpasssummary(5) = test_pass;

%% Step 4 - Cleanup
clearvars -except  curDir  ConfigDir ...
    ntests testpasssummary...
    plot2vtk VTKFolder progress_on...
    saveGraphs_on saveMatData_on plotGraphs_on saveVideo_on