% ------------------------------------------------------------------------
% Runs unit Test 11 - L2/L3 manufactured solution convergence as a part of RunTests
% ------------------------------------------------------------------------
% Calculates the convergence rates of a uniform L3/L2 mesh using a
% manufactured solution in which
% u = sin(xt)
% p = sin(xt)
% Both fields are prescribed, but the problem is uncoupled
% ------------------------------------------------------------------------
% Adapted from https://github.com/GCMLab (Acknowledgements: Bruce Gee)
% ------------------------------------------------------------------------

% Create test VTK folder
if plot2vtk
    vtk_dir = fullfile(VTKFolder,'\Test11');
    if ~isfolder(vtk_dir)
        mkdir(vtk_dir)
    end
end
% test runs 3 meshes, only finest mesh will be saved

fprintf('\n\n Test 11: Manufactured Solution - transient, L3/L2 elements\n')

%% Step 1 - Run Simulation
config_name = 'ManufacturedSolutionTransient';
meshfilename = '';

% Run coarse mesh
nelements = 16;
main
% store variables coarse mesh
u_coarse = Solution.u;
stress_coarse = stress;
strain_coarse = strain;
p_coarse = Solution.p;
gradp_coarse = gradp;
flux_coarse = flux;
MeshU_coarse = MeshU;
MeshP_coarse = MeshP;

% Run fine mesh
nelements = 32;
main
% store variables fine mesh
u_fine = Solution.u;
stress_fine = stress;
strain_fine = strain;
p_fine = Solution.p;
gradp_fine = gradp;
flux_fine = flux;
MeshU_fine = MeshU;
MeshP_fine = MeshP;

% Run finer mesh
nelements = 64;
main
% store variables finer mesh
u_finer = Solution.u;
stress_finer = stress;
strain_finer = strain;
p_finer = Solution.p;
gradp_finer = gradp;
flux_finer = flux;
MeshU_finer = MeshU;
MeshP_finer = MeshP;

%% Step 2 - Check results
% u field
BC.ux = @(x) sin(x*Control.tend);
BC.dudx = @(x) cos(x*Control.tend);
[m_L2u, m_eu] = ManufacturedSolution1D_check(u_coarse, u_fine, u_finer, ...
    stress_coarse, stress_fine, stress_finer, strain_coarse, strain_fine, ...
    strain_finer, MeshU_coarse, MeshU_fine, MeshU_finer, Material, Control, QuadU, BC);
% p field
BC.ux = @(x) sin(x*Control.tend);
BC.dudx = @(x) cos(x*Control.tend);
[m_L2p, m_ep] = ManufacturedSolution1D_check(p_coarse, p_fine, p_finer, ...
    flux_coarse, flux_fine, flux_finer, gradp_coarse, gradp_fine, ...
    gradp_finer, MeshP_coarse, MeshP_fine, MeshP_finer, Material, Control, QuadU, BC);

fprintf('\nL3 L2-norm converges at a rate of %.2f for u',m_L2u);
fprintf('\nL3  e-norm converges at a rate of %.2f for u',m_eu);
fprintf('\nL3 L2-norm converges at a rate of %.2f for p',m_L2p);
fprintf('\nL3  e-norm converges at a rate of %.2f for p',m_ep);

convergence_tolerance = 0.05;
if m_L2u >= (3 - convergence_tolerance) && m_eu >= (2 - convergence_tolerance) && m_L2p >= (2 - convergence_tolerance) && m_ep >= (1 - convergence_tolerance)
    test_pass = 1;
else
    test_pass = 0;
end

%% Step 3 - Output results
if test_pass
    fprintf('\nPASS Convergence Test Transient\n')
else
    fprintf('\n\nFAIL Convergence Test Transient\n')
end
testpasssummary(10) = test_pass;

%% Step 4 - Cleanup
clearvars -except  curDir  ConfigDir ...
    ntests testpasssummary...
    plot2vtk VTKFolder progress_on...
    saveGraphs_on saveMatData_on plotGraphs_on