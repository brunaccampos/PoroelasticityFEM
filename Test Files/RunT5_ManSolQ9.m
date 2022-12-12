% ------------------------------------------------------------------------
% Runs unit Test 5 - Q9 manufactured solution convergence as a part of RunTests
% ------------------------------------------------------------------------
% Test 5 calculates the convergence rates of a uniform Q9 mesh using a
% manufactured solution in which
% ux = x^5 + x*y^3 - y^6
% uy = x^5 + x*y^3 - y^6
% under plane stress conditions
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
% test runs 3 meshes, only finest mesh will be saved

fprintf('\n\n Test 5: Manufactured Solution - Q9 elements\n')

%% Step 1 - Run Simulation
config_name = 'ManufacturedSolutionQ9';

% Run coarse mesh
meshfilename = 'Mesh Files\Manufactured_coarseQ9.msh';
mainTests
% store variables coarse mesh
d_coarse = u;
stress_coarse = stress;
strain_coarse = strain;
Mesh_coarse = MeshU;

% Run fine mesh
meshfilename = 'Mesh Files\Manufactured_fineQ9.msh';
mainTests
% store variables fine mesh
d_fine = u;
stress_fine = stress;
strain_fine = strain;
Mesh_fine = MeshU;

% Run finer mesh
meshfilename = 'Mesh Files\Manufactured_finerQ9.msh';
mainTests
% store variables finer mesh
d_finer = u;
stress_finer = stress;
strain_finer = strain;
Mesh_finer = MeshU;

%% Step 2 - Check results
[m_L2, m_e] = ManufacturedSolution_check(d_coarse, d_fine, d_finer, ...
    stress_coarse, stress_fine, stress_finer, strain_coarse, strain_fine, ...
    strain_finer, Mesh_coarse, Mesh_fine, Mesh_finer, Material, Control, Quad);

fprintf('\nQ9 L2-norm converges at a rate of %.2f',m_L2)
fprintf('\nQ9  e-norm converges at a rate of %.2f',m_e)

convergence_tolerance = 0.05;
if m_L2 >= (3 - convergence_tolerance) && m_e >= (2 - convergence_tolerance)
    test_pass = 1;
else
    test_pass = 0;
end

%% Step 3 - Output results
if test_pass
    fprintf('\nPASS\n')
else
    fprintf('\n\nFAIL\n')
end
testpasssummary(5) = test_pass;

%% Step 4 - Cleanup
clearvars -except  curDir  ConfigDir ...
    ntests testpasssummary...
    plot2vtk VTKFolder progress_on