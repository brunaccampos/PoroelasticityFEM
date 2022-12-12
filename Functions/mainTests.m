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
[Material, MeshU, MeshP, BC, Control] = feval(config_name, ConfigDir, progress_on, meshfilename);

%% Quadrature points
Quad = GlobalQuad(MeshU, Control);

%% Assemble system matrices
disp([num2str(toc),': Assembling System Matrices...']);
[Kuu, Kup, Kpp, S] = ComputeSystemMatrices_Steady(Material, MeshU, MeshP, Control, Quad);

%% Assemble system load vectors
[fu,fp] = ComputeSystemLoads(BC, MeshU, MeshP, Control, Quad);

%% Initialize iteration variables
nsd = MeshU.nsd; % number of spatial directions

%% Solve system
[u,~, p, ~,fE, ~] = SolverSteady(Kuu, Kup, Kpp, S, fu, fp, BC, Control, []);
fu(BC.fixed_u) = fE;

%% Compute strains and stresses
[strain, stress] = ComputeSolidStress(Material, MeshU, u);