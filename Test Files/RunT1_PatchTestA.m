% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later
% ------------------------------------------------------------------------
% Runs unit Test 1 - Patch Test A as part of RunTests
% ------------------------------------------------------------------------
% For Patch Test A, all nodes are restrained and nodal displacement values
% are specfied according to the exact solution. The error between the FEA
% and exact solutions is then calculated. The FEA approximate solution
% should be exact.
% ------------------------------------------------------------------------

% Create test VTK folder
if plot2vtk
    vtk_dir = fullfile(VTKFolder,'\Test1');
    if ~isfolder(vtk_dir)
        mkdir(vtk_dir)
    end
end

fprintf('\n\n Test 1: Patch Test A - Q4 elements\n')

%% Step 1 - Run Simulation
config_name = 'PatchTestA';
meshfilename = 'Mesh Files\PatchTest.msh';
nelements = 0;
main

%% Step 2 - Check results
[disp_er, stress_er, reaction_er] = PatchTest_check(Solution.u, stress, fu, MeshU, BC, Material);

fprintf('\nQ4-patch test A: Displacement error is %.2f',disp_er)
fprintf('\nQ4-patch test A: Stress error is %.2f',stress_er)
fprintf('\nQ4-patch test A: Reaction forces error is %.2f',reaction_er)

convergence_tolerance = 1e-10;
if disp_er <= convergence_tolerance && stress_er <= convergence_tolerance && reaction_er <= convergence_tolerance
    test_pass = 1;
else
    test_pass = 0;
end

%% Step 3 - Output results
if test_pass
    fprintf('\nPASS Patch Test A\n')
else
    fprintf('\n\nFAIL Patch Test A\n')
end
testpasssummary(1) = test_pass;

%% Step 4 - Cleanup
clearvars -except  curDir  ConfigDir ...
    ntests testpasssummary...
    plot2vtk VTKFolder progress_on...
    saveGraphs_on saveMatData_on plotGraphs_on saveVideo_on