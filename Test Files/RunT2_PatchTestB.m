% ------------------------------------------------------------------------
% Runs unit Test 2 - Patch Test B as part of RunTests
% ------------------------------------------------------------------------
% For Patch Test B, only nodes 1-8 (nodes in the boundaries) are restrained
% with their displacements specified according to the exact solution. The error between the FEA
% and exact solutions is then calculated. The FEA approximate solution
% should be exact.

% Adapted from https://github.com/GCMLab (Acknowledgements: Bruce Gee)
% ------------------------------------------------------------------------

% Create test VTK folder
if plot2vtk
    vtk_dir = fullfile(VTKFolder,'\Test2');
    if ~isfolder(vtk_dir)
        mkdir(vtk_dir)
    end
end

fprintf('\n\n Test 2: Patch Test B - Q4 elements\n')

%% Step 1 - Run Simulation
config_name = 'PatchTestB';
meshfilename = 'Mesh Files\PatchTest_T.msh';
main

%% Step 2 - Check results
[disp_er, stress_er, reaction_er] = PatchTest_check(Solution.u, stress, fu, MeshU, BC, Material);

fprintf('\nQ4-patch test B: Displacement error is %.2f',disp_er)
fprintf('\nQ4-patch test B: Stress error is %.2f',stress_er)
fprintf('\nQ4-patch test B: Reaction forces error is %.2f',reaction_er)

convergence_tolerance = 1e-10;
if disp_er <= convergence_tolerance && stress_er <= convergence_tolerance && reaction_er <= convergence_tolerance
    test_pass = 1;
else
    test_pass = 0;
end

%% Step 3 - Output results
if test_pass
    fprintf('\nPASS Patch Test B\n')
else
    fprintf('\n\nFAIL Patch Test B\n')
    return
end
testpasssummary(2) = test_pass;

%% Step 4 - Cleanup
clearvars -except  curDir  ConfigDir ...
    ntests testpasssummary...
    plot2vtk VTKFolder progress_on...
    saveGraphs_on saveMatData_on plotGraphs_on