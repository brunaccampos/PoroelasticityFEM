function PostProcessingModes(mode, phi_u, phi_p, MeshU, MeshP, BC, config_name, vtk_dir)
% ------------------------------------------------------------------------
% Export results to VTK file
% ------------------------------------------------------------------------
% Adapted from: https://github.com/GCMLab (Acknowledgements: Matin Parchei Esfahani)

%% Initialize variables
nsd = MeshU.nsd; % number of spatial dimensions
plotmodeu = zeros(MeshU.nDOF,1);
plotmodeu(BC.free_u) = phi_u(:, mode);
plotmodep = zeros(MeshP.nDOF,1);
plotmodep(BC.free_p) = phi_p(:, mode);

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
    % mode shape displacement
    scalardataU(1).name = 'mode_u';
    scalardataU(1).data = plotmodeu(xdofs_u);
    scalardataU(1).type = 'float';
    % fixed DOFs
    scalardataU(2).name = 'fixedU';  scalardataU(2).data = fixedU(MeshU.DOF);
    scalardataU(2).type = 'int';

    vectordataU = [];
    
elseif nsd == 2
    xdofs_u = MeshU.DOF(:,1); % DOFs in x
    ydofs_u = MeshU.DOF(:,2); % DOFs in y
    
    %% Storing data for solid media - 2D transient case  
    % mode shape displacement
    vectordataU(1).name = 'mode_u';
    vectordataU(1).data = [plotmodeu(xdofs_u) plotmodeu(ydofs_u) zeros(length(xdofs_u),1)];
    vectordataU(1).type = 'float';
    
    % fixed DOFs
    scalardataU(1).name = 'fixedU';  scalardataU(1).data = fixedU(MeshU.DOF);
    scalardataU(1).type = 'int';
end

%% Storing data for fluid media
% pressure
scalardataP(1).name = 'mode_p';
scalardataP(1).data = plotmodep;
scalardataP(1).type = 'float';

vectordataP = [];

%% Write to VTK
% solid
description = config_name; % config file name
nameU = 'Mode_U.vtk.';
filenameU = fullfile(vtk_dir, [nameU, num2str(mode)]);
WriteMesh2VTK(filenameU, description, MeshU, scalardataU, vectordataU);

% fluid
description = config_name; % config file name
nameP = 'Mode_P.vtk.';
filenameP = fullfile(vtk_dir, [nameP, num2str(mode)]);
WriteMesh2VTK(filenameP, description, MeshP, scalardataP, vectordataP);


end