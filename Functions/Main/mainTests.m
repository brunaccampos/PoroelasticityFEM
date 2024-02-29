% Porous Media Simulation - Tests
% ------------------------------------------------------------------------
% Created by Bruna Campos
% bccampos@uwaterloo.ca
% Department of Civil Engineering, University of Waterloo
% January 2022
% ------------------------------------------------------------------------
% Reference: https://github.com/GCMLab
% ------------------------------------------------------------------------

%% Delete past files
if plot2vtk
    if progress_on
        fprintf('%.2f: Deleting old vtk and csv files...\n', toc);
    end
    sol = fullfile(vtk_dir, '*');
    delete(sol)
end

%% Input problem
if progress_on
    fprintf('%.2f: Reading config file...\n', toc);
end
[Material, MeshU, MeshP, MeshN, BC, Control] = feval(config_name, ConfigDir, progress_on, meshfilename);

%% Set defaults
[Material, BC, Control] = setDefaults(Material, MeshU, MeshP, BC, Control);

%% Quadrature points
QuadU = GlobalQuad(MeshU, Control);
QuadP = GlobalQuad(MeshP, Control);

%% Assemble system matrices
disp([num2str(toc),': Assembling System Matrices...']);
[Kuu, Kup, Kpp, Kpu, S] = ComputeMatricesTr1_Biot_UP(Material, MeshU, MeshP, QuadU, QuadP);

%% Assemble system load vectors
[fu,fp,~] = ComputeLoads(BC, MeshU, MeshP, [], Control, QuadU, QuadP);

%% Solve system
Solution = SolverTr_UP(Kuu, Kup, Kpp, Kpu, S, fu, fp, BC, Control, []);
fu(BC.fixed_u) = Solution.fE;
fp(BC.fixed_p) = Solution.qE;

%% Compute strains and stresses
[strain, stress] = ComputeSolidStress(Material, MeshU, Solution.u);
[gradp, flux] = ComputeFluidFlux(Material, MeshP, Solution.p);

% post processing: export VTK file
if plot2vtk
    PostProcessing(Solution, Material, MeshU, MeshP, [], Control, BC, config_name, vtk_dir);
end
