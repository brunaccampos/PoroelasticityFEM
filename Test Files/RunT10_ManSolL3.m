% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later
% ------------------------------------------------------------------------
% Runs unit Test 10 - L3 manufactured solution convergence as a part of RunTests
% ------------------------------------------------------------------------
% Calculates the convergence rates of a uniform L3 mesh using a
% manufactured solution in which
% ux = x^5 - x^4
% ------------------------------------------------------------------------

% Create test VTK folder
if plot2vtk
    vtk_dir = fullfile(VTKFolder,'\Test10');
    if ~isfolder(vtk_dir)
        mkdir(vtk_dir)
    end
end
% test runs 3 meshes, only finest mesh will be saved

fprintf('\n\n Test 10: Manufactured Solution - L3 elements\n')

%% Step 1 - Run Simulation
config_name = 'ManufacturedSolutionL3';
meshfilename = '';

% Run coarse mesh
nelements = 16;
main
% store variables coarse mesh
d_coarse = Solution.u;
stress_coarse = stress;
strain_coarse = strain;
Mesh_coarse = MeshU;

% Run fine mesh
nelements = 32;
main
% store variables fine mesh
d_fine = Solution.u;
stress_fine = stress;
strain_fine = strain;
Mesh_fine = MeshU;

% Run finer mesh
nelements = 64;
main
% store variables finer mesh
d_finer = Solution.u;
stress_finer = stress;
strain_finer = strain;
Mesh_finer = MeshU;

%% Step 2 - Check results
[m_L2, m_e] = ManufacturedSolution1D_check(d_coarse, d_fine, d_finer, ...
    stress_coarse, stress_fine, stress_finer, strain_coarse, strain_fine, ...
    strain_finer, Mesh_coarse, Mesh_fine, Mesh_finer, Material, Control, QuadU, BC);

fprintf('\nL3 L2-norm converges at a rate of %.2f',m_L2)
fprintf('\nL3  e-norm converges at a rate of %.2f',m_e)

convergence_tolerance = 0.05;
if m_L2 >= (3 - convergence_tolerance) && m_e >= (2 - convergence_tolerance)
    test_pass = 1;
else
    test_pass = 0;
end

%% Step 3 - Output results
if test_pass
    fprintf('\nPASS Convergence Test L3\n')
else
    fprintf('\n\nFAIL Convergence Test L3\n')
end
testpasssummary(10) = test_pass;

%% Step 4 - Cleanup
clearvars -except  curDir  ConfigDir ...
    ntests testpasssummary...
    plot2vtk VTKFolder progress_on...
    saveGraphs_on saveMatData_on plotGraphs_on saveVideo_on