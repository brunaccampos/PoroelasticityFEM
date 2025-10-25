% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later
% ------------------------------------------------------------------------
% Runs unit Test 3 - Patch Test C as part of RunTests
% ------------------------------------------------------------------------
% Patch Test C is performed with node 1 fully restrained and nodes 4 and 8
% restrained only in the x -direction. Nodal forces are applied to nodes 2,
% 3, and 6 in accordance with the values generated through the boundary
% tractions by sigma(x)=2. The error between the FEA
% and exact solutions is then calculated. The FEA approximate solution
% should be exact.
% ------------------------------------------------------------------------

% Create test VTK folder
if plot2vtk
    vtk_dir = fullfile(VTKFolder,'\Test3');
    if ~isfolder(vtk_dir)
        mkdir(vtk_dir)
    end
end

fprintf('\n\n Test 3: Patch Test C - Q4 elements\n')

%% Step 1 - Run Simulation
config_name = 'PatchTestC';
meshfilename = 'Mesh Files\PatchTest.msh';
nelements = 0;
main

%% Step 2 - Check results
[disp_er, stress_er, reaction_er] = PatchTest_check(Solution.u, stress, fu, MeshU, BC, Material);

fprintf('\nQ4-patch test C: Displacement error is %.2f',disp_er)
fprintf('\nQ4-patch test C: Stress error is %.2f',stress_er)
fprintf('\nQ4-patch test C: Reaction forces error is %.2f',reaction_er)

convergence_tolerance = 1e-10;
if disp_er <= convergence_tolerance && stress_er <= convergence_tolerance && reaction_er <= convergence_tolerance
    test_pass = 1;
else
    test_pass = 0;
end

%% Step 3 - Output results
if test_pass
    fprintf('\nPASS Patch Test C\n')
else
    fprintf('\n\nFAIL Patch Test C\n')
    return
end
testpasssummary(3) = test_pass;

%% Step 4 - Cleanup
clearvars -except  curDir  ConfigDir ...
    ntests testpasssummary...
    plot2vtk VTKFolder progress_on...
    saveGraphs_on saveMatData_on plotGraphs_on saveVideo_on