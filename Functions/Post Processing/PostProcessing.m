function PostProcessing(Solution, Material, MeshU, MeshP, MeshN, Control, BC, config_name, vtk_dir)
% ------------------------------------------------------------------------
% Export results to VTK file
% ------------------------------------------------------------------------
% Adapted from: https://github.com/GCMLab (Acknowledgements: Matin Parchei Esfahani)

%% Initialize variables
step = Control.step; % time step number
% step = 0; % for running patch tests

nsd = MeshU.nsd; % number of spatial dimensions

%% Compute solid stress
[strainU, stressU] = ComputeSolidStress(Material, MeshU, Solution.u);

%% Compute fluid flux
[q] = ComputeFluidFlux(Material, MeshP, Solution.p);

%% Store fixed DOFs
% fixed DOFs displacement
fixedU = zeros(MeshU.nDOF,1);
fixedU(BC.fixed_u) = 1;
% fixed DOFs pressure
fixedP = zeros(MeshP.nDOF,1);
fixedP(BC.fixed_p) = 1;

if nsd == 1
    xdofs_u = MeshU.DOF(:,1); % DOFs
    
    %% Storing data for solid media - 1D case
    % displacement
    scalardataU(1).name = 'disp_u';
    scalardataU(1).data = Solution.u(xdofs_u);
    scalardataU(1).type = 'float';
    % velocity
    scalardataU(2).name = 'vel_u';
    scalardataU(2).data = Solution.udot(xdofs_u);
    scalardataU(2).type = 'float';
    % strain
    scalardataU(3).name = 'exx';     scalardataU(3).data = strainU(:,1);
    scalardataU(3).type = 'float';
    % stress
    scalardataU(4).name = 'Sxx';     scalardataU(4).data = stressU(:,1);
    scalardataU(4).type = 'float';
    % fixed DOFs
    scalardataU(5).name = 'fixedU';  scalardataU(5).data = fixedU(MeshU.DOF);
    scalardataU(5).type = 'int';
    if contains(Control.Biotmodel, 'Dynamic')
        % acceleration
        scalardataU(6).name = 'acc_u';
        scalardataU(6).data = Solution.u2dot(xdofs_u);
        scalardataU(6).type = 'float';
    end
    vectordataU = [];
    
elseif nsd == 2
    xdofs_u = MeshU.DOF(:,1); % DOFs in x
    ydofs_u = MeshU.DOF(:,2); % DOFs in y
    
    %% Storing data for solid media - 2D transient case  
    % displacement
    vectordataU(1).name = 'disp_u';
    vectordataU(1).data = [Solution.u(xdofs_u) Solution.u(ydofs_u) zeros(length(xdofs_u),1)];
    vectordataU(1).type = 'float';
    % velocity
    vectordataU(2).name = 'vel_u';
    vectordataU(2).data = [Solution.udot(xdofs_u) Solution.udot(ydofs_u) zeros(length(xdofs_u),1)];
    vectordataU(2).type = 'float';
    if contains(Control.Biotmodel, 'Dynamic')
        % acceleration
        vectordataU(3).name = 'acc_u';
        vectordataU(3).data = [Solution.u2dot(xdofs_u) Solution.u2dot(ydofs_u) zeros(length(xdofs_u),1)];
        vectordataU(3).type = 'float';
    end
    % strain
    scalardataU(1).name = 'exx';     scalardataU(1).data = strainU(:,1);
    scalardataU(1).type = 'float';
    scalardataU(2).name = 'eyy';     scalardataU(2).data = strainU(:,2);
    scalardataU(2).type = 'float';
    scalardataU(3).name = 'exy';     scalardataU(3).data = stressU(:,3);
    scalardataU(3).type = 'float';
    % stress
    scalardataU(4).name = 'Sxx';     scalardataU(4).data = stressU(:,1);
    scalardataU(4).type = 'float';
    scalardataU(5).name = 'Syy';     scalardataU(5).data = stressU(:,2);
    scalardataU(5).type = 'float';
    scalardataU(6).name = 'Sxy';     scalardataU(6).data = stressU(:,3);
    scalardataU(6).type = 'float';
    % fixed DOFs
    scalardataU(7).name = 'fixedU';  scalardataU(7).data = fixedU(MeshU.DOF);
    scalardataU(7).type = 'int';
    
end

%% Storing data for fluid media
% pressure
scalardataP(1).name = 'pressure';
scalardataP(1).data = Solution.p;
scalardataP(1).type = 'float';
% fixed DOFs 
scalardataP(2).name = 'fixedP';  scalardataP(2).data = fixedP(MeshP.DOF);
scalardataP(2).type = 'int';

if nsd == 1
    % fluid flux in 1D
    scalardataP(4).name = 'q';
    scalardataP(4).data = q;
    scalardataP(4).type = 'float';
    vectordataP = [];
elseif nsd == 2
    % fluid flux in 2D
    vectordataP(1).name = 'q';
    vectordataP(1).data = q;
    vectordataP(1).type = 'float';
end

%% Storing data for porosity field
% porosity
if contains(Control.Biotmodel, 'Spanos')
    scalardataN(1).name = 'porosity';
    scalardataN(1).data = Solution.n;
    scalardataN(1).type = 'float';
else
    scalardataN = [];
end

%% Write to VTK
% solid
description = config_name; % config file name
nameU = 'Solution_U.vtk.';
filenameU = fullfile(vtk_dir, [nameU, num2str(step)]);
WriteMesh2VTK(filenameU, description, MeshU, scalardataU, vectordataU);

% fluid
description = config_name; % config file name
nameP = 'Solution_P.vtk.';
filenameP = fullfile(vtk_dir, [nameP, num2str(step)]);
WriteMesh2VTK(filenameP, description, MeshP, scalardataP, vectordataP);

% porosity
if contains(Control.Biotmodel, 'Spanos')
    description = config_name; % config file name
    nameN = 'Solution_N.vtk.';
    filenameN = fullfile(vtk_dir, [nameN, num2str(step)]);
    WriteMesh2VTK(filenameN, description, MeshN, scalardataN, []);
end

end