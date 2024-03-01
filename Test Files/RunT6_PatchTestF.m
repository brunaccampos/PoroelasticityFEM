% ------------------------------------------------------------------------
% Runs unit Test 6 - Patch Test F as part of RunTests
% ------------------------------------------------------------------------
% Patch Test C is performed with node 1 fully restrained and nodes 4 and 8
% restrained only in the x -direction. Nodal fluxes are applied to nodes 2,
% 3, and 6 in accordance with the values generated through the boundary
% fluxes. The error between the FEA and exact solutions is then calculated. 
% The FEA approximate solution should be exact.
% ------------------------------------------------------------------------
% Adapted from https://github.com/GCMLab (Acknowledgements: Bruce Gee)
% ------------------------------------------------------------------------

% Create test VTK folder
if plot2vtk
    vtk_dir = fullfile(VTKFolder,'\Test6');
    if ~isfolder(vtk_dir)
        mkdir(vtk_dir)
    end
end

fprintf('\n\n Test 6: Patch Test F - Q4 elements\n')

%% Step 1 - Run Simulation
config_name = 'PatchTestF';
meshfilename = 'Mesh Files\PatchTest.msh';
main

%% Step 2 - Check results
[press_er, flux_er] = PatchTest_check_v2(Solution.p, flux, MeshP, BC, Material);

fprintf('\nQ4-patch test F: Pressure error is %.2f',press_er)
fprintf('\nQ4-patch test F: Flux error is %.2f',flux_er)

convergence_tolerance = 5e-10;
if press_er <= convergence_tolerance && flux_er <= convergence_tolerance
    test_pass = 1;
else
    test_pass = 0;
end

%% Step 3 - Output results
if test_pass
    fprintf('\nPASS Patch Test F\n')
else
    fprintf('\n\nFAIL Patch Test F\n')
    return
end
testpasssummary(6) = test_pass;

%% Step 4 - Cleanup
clearvars -except  curDir  ConfigDir ...
    ntests testpasssummary...
    plot2vtk VTKFolder progress_on...
    saveGraphs_on saveMatData_on plotGraphs_on