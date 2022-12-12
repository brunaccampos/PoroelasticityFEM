function PostProcessing(u, udot, u2dot, p, Material, MeshU, MeshP, Control, Quad, BC, config_name, vtk_dir)
% Export results to VTK file
% ------------------------------------------------------------------------
% Input
% ------------------------------------------------------------------------

% ------------------------------------------------------------------------
% Adapted from: https://github.com/GCMLab (Acknowledgements: Matin Parchei Esfahani)

%% Initialize variables
step = Control.step; % time step number
nsd = MeshU.nsd; % number of spatial dimensions

%% Compute solid stress
[strainU, stressU] = ComputeSolidStress(Material, MeshU, Control, Quad, u);

% fixed DOFs displacement
fixedU = zeros(MeshU.nDOF,1);
fixedU(BC.fixed_u) = 1;
% fixed DOFs pressure
fixedP = zeros(MeshP.nDOF,1);
fixedP(BC.fixed_p) = 1;
% fixed DOFs pressure gradient
fixedfP = zeros(MeshP.nDOF,1);
fixedfP(BC.fixed_fp) = 1;

if nsd == 1
    xdofs_u = MeshU.DOF(:,1); % DOFs

    %% Storing data for solid media - 1D case
    % displacement
    scalardataU(1).name = 'disp_u';
    scalardataU(1).data = u(xdofs_u);
    scalardataU(1).type = 'float';
    % velocity
    scalardataU(2).name = 'vel_u';
    scalardataU(2).data = udot(xdofs_u);
    scalardataU(2).type = 'float';
    % acceleration
    scalardataU(3).name = 'acc_u';
    scalardataU(3).data = u2dot(xdofs_u);
    scalardataU(3).type = 'float';
    % strain
    scalardataU(4).name = 'exx';     scalardataU(4).data = strainU(:,1);
    scalardataU(4).type = 'float';
    % stress
    scalardataU(5).name = 'Sxx';     scalardataU(5).data = stressU(:,1);
    scalardataU(5).type = 'float';
    % fixed DOFs
    scalardataU(6).name = 'fixedU';  scalardataU(6).data = fixedU(MeshU.DOF);
    scalardataU(6).type = 'int';

elseif nsd == 2
    xdofs_u = MeshU.DOF(:,1); % DOFs in x
    ydofs_u = MeshU.DOF(:,2); % DOFs in y

    %% Storing data for solid media - 2D case
    % displacement
    scalardataU(1).name = 'disp_u';
    scalardataU(1).data = [u(xdofs_u) u(ydofs_u) zeros(length(xdofs_u),1)];
    scalardataU(1).type = 'float';
    % velocity
    scalardataU(2).name = 'vel_u';
    scalardataU(2).data = [udot(xdofs_u) udot(ydofs_u) zeros(length(xdofs_u),1)];
    scalardataU(2).type = 'float';
    % acceleration
    scalardataU(3).name = 'acc_u';
    scalardataU(3).data = [u2dot(xdofs_u) u2dot(ydofs_u) zeros(length(xdofs_u),1)];
    scalardataU(3).type = 'float';
    % strain
    scalardataU(4).name = 'exx';     scalardataU(4).data = strainU(:,1);
    scalardataU(4).type = 'float';
    scalardataU(5).name = 'eyy';     scalardataU(5).data = strainU(:,2);
    scalardataU(5).type = 'float';
    scalardataU(6).name = 'exy';     scalardataU(6).data = stressU(:,3);
    scalardataU(6).type = 'float';
    % stress
    scalardataU(7).name = 'Sxx';     scalardataU(7).data = stressU(:,1);
    scalardataU(7).type = 'float';
    scalardataU(8).name = 'Syy';     scalardataU(8).data = stressU(:,2);
    scalardataU(8).type = 'float';
    scalardataU(9).name = 'Sxy';     scalardataU(9).data = stressU(:,3);
    scalardataU(9).type = 'float';
    % fixed DOFs
    scalardataU(10).name = 'fixedU';  scalardataU(10).data = fixedU(MeshU.DOF);
    scalardataU(10).type = 'int';
end

%% Storing data for fluid media
% pressure
scalardataP(1).name = 'pressure';
scalardataP(1).data = p;
scalardataP(1).type = 'float';

% fixed DOFs
scalardataP(2).name = 'fixedP';  scalardataP(2).data = fixedP(MeshP.DOF);
scalardataP(2).type = 'int';
% fixed pressure
scalardataP(3).name = 'fixedFP';  scalardataP(3).data = fixedfP(MeshP.DOF);
scalardataP(3).type = 'int';

%% Write to VTK
% solid
description = config_name; % config file name
nameU = 'Solution_U.vtk.';
filenameU = fullfile(vtk_dir, [nameU, num2str(step)]);

WriteMesh2VTK(filenameU, description, MeshU, scalardataU);

% fluid
description = config_name; % config file name
nameP = 'Solution_P.vtk.';
filenameP = fullfile(vtk_dir, [nameP, num2str(step)]);

WriteMesh2VTK(filenameP, description, MeshP, scalardataP);

end