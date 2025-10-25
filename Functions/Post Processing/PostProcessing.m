% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function PostProcessing(Solution, Material, MeshU, MeshP, MeshN, Control, BC, config_name, vtk_dir)
% Export results to VTK file
% Acknowledgements: Matin Parchei Esfahani

%% Initialize variables
step = Control.step; % time step number
nsd = MeshU.nsd; % number of spatial dimensions

%% Compute solid stress
[strainU, stressU] = ComputeSolidStress(Material, MeshU, Solution.u);

%% Compute fluid flux
[q] = ComputeFluidFlux(Material, MeshP, Solution.p);

%% Compute porosity
if contains(Control.PMmodel, 'UPU') || contains(Control.PMmodel, 'UPV') || contains(Control.PMmodel, 'UPW')
    [eta, etadot] = ComputePorosity(Material, MeshU, Solution, Control);
end

%% Store fixed DOFs
% fixed DOFs solid displacement
fixedU = zeros(MeshU.nDOF,1);
fixedU(BC.fixed_u) = 1;
% fixed DOFs fluid displacement
if contains(Control.PMmodel, 'UPU')
    fixedUf = zeros(MeshU.nDOF,1);
    fixedUf(BC.fixed_uf) = 1;
end
% fixed DOFs fluid velocity
if contains(Control.PMmodel, 'UPV')
    fixedUfdot = zeros(MeshU.nDOF,1);
    fixedUfdot(BC.fixed_ufdot) = 1;
end
% fixed DOFs relative fluid velocity
if contains(Control.PMmodel, 'UPW')
    fixedW = zeros(MeshU.nDOF,1);
    fixedW(BC.fixed_w) = 1;
end
% fixed DOFs pressure
fixedP = zeros(MeshP.nDOF,1);
fixedP(BC.fixed_p) = 1;

if nsd == 1
    xdofs_u = MeshU.DOF(:,1); % DOFs
    
    %% Solid field - 1D
    % fixed DOFs
    scalardataU(1).name = 'fixed_us';
    scalardataU(1).data = fixedU(MeshU.DOF);
    scalardataU(1).type = 'int';
    % displacement
    scalardataU(end+1).name = 'disp_us';
    scalardataU(end).data = Solution.u(xdofs_u);
    scalardataU(end).type = 'float';
    % velocity
    scalardataU(end+1).name = 'vel_us';
    scalardataU(end).data = Solution.udot(xdofs_u);
    scalardataU(end).type = 'float';
       
    if contains(Control.PMmodel, 'Dyn')
        % acceleration
        scalardataU(end+1).name = 'acc_us';
        scalardataU(end).data = Solution.u2dot(xdofs_u);
        scalardataU(end).type = 'float';
    end
    
    % strain
    scalardataU(end+1).name = 'exx_us';
    scalardataU(end).data = strainU(:,1);
    scalardataU(end).type = 'float';
    % stress
    scalardataU(end+1).name = 'Sxx_us';
    scalardataU(end).data = stressU(:,1);
    scalardataU(end).type = 'float';
    
    vectordataU = [];
    
    %% Fluid displacement - 1D
    if contains(Control.PMmodel, 'UPU')
        % fixed DOFs
        scalardataUf(1).name = 'fixed_uf';
        scalardataUf(1).data = fixedUf(MeshU.DOF);
        scalardataUf(1).type = 'int';
        % displacement
        scalardataUf(end+1).name = 'disp_uf';
        scalardataUf(end).data = Solution.uf(xdofs_u);
        scalardataUf(end).type = 'float';
        % velocity
        scalardataUf(end+1).name = 'vel_uf';
        scalardataUf(end).data = Solution.ufdot(xdofs_u);
        scalardataUf(end).type = 'float';
        
        if contains(Control.PMmodel, 'Dyn')
            % acceleration
            scalardataUf(end+1).name = 'acc_uf';
            scalardataUf(end).data = Solution.uf2dot(xdofs_u);
            scalardataUf(end).type = 'float';
        end
        vectordataUf = [];
    end
    
    %% Fluid velocity - 1D
    if contains(Control.PMmodel, 'UPV')
        % fixed DOFs
        scalardataUfdot(1).name = 'fixed_ufdot';
        scalardataUfdot(1).data = fixedUfdot(MeshU.DOF);
        scalardataUfdot(1).type = 'int';
        % velocity
        scalardataUfdot(end+1).name = 'vel_uf';
        scalardataUfdot(end).data = Solution.ufdot(xdofs_u);
        scalardataUfdot(end).type = 'float';
        
        if contains(Control.PMmodel, 'Dyn')
            % acceleration
            scalardataUfdot(end+1).name = 'acc_uf';
            scalardataUfdot(end).data = Solution.uf2dot(xdofs_u);
            scalardataUfdot(end).type = 'float';
        end
        vectordataUfdot = [];
    end
    
    %% Fluid velocity (relative) - 1D
    if contains(Control.PMmodel, 'UPW')
        % fixed DOFs
        scalardataW(1).name = 'fixed_w';
        scalardataW(1).data = fixedW(MeshU.DOF);
        scalardataW(1).type = 'int';
        % velocity
        scalardataW(end+1).name = 'vel_w';
        scalardataW(end).data = Solution.w(xdofs_u);
        scalardataW(end).type = 'float';
        
        if contains(Control.PMmodel, 'Dyn')
            % acceleration
            scalardataW(end+1).name = 'acc_w';
            scalardataW(end).data = Solution.wdot(xdofs_u);
            scalardataW(end).type = 'float';
        end
        vectordataW = [];
    end

elseif nsd == 2
    xdofs_u = MeshU.DOF(:,1); % DOFs in x
    ydofs_u = MeshU.DOF(:,2); % DOFs in y
    
    %% Storing data for solid media - 2D case
    % fixed DOFs
    scalardataU(1).name = 'fixed_us';
    scalardataU(1).data = fixedU(MeshU.DOF);
    scalardataU(1).type = 'int';
    % displacement
    vectordataU(1).name = 'disp_us';
    vectordataU(1).data = [Solution.u(xdofs_u) Solution.u(ydofs_u) zeros(length(xdofs_u),1)];
    vectordataU(1).type = 'float';
    % velocity
    vectordataU(end+1).name = 'vel_us';
    vectordataU(end).data = [Solution.udot(xdofs_u) Solution.udot(ydofs_u) zeros(length(xdofs_u),1)];
    vectordataU(end).type = 'float';
    
    if contains(Control.PMmodel, 'Dyn')
        % acceleration
        vectordataU(end+1).name = 'acc_us';
        vectordataU(end).data = [Solution.u2dot(xdofs_u) Solution.u2dot(ydofs_u) zeros(length(xdofs_u),1)];
        vectordataU(end).type = 'float';
    end
    
    % strain xx
    scalardataU(end+1).name = 'exx_us';
    scalardataU(end).data = strainU(:,1);
    scalardataU(end).type = 'float';
    % strain yy
    scalardataU(end+1).name = 'eyy_us';
    scalardataU(end).data = strainU(:,2);
    scalardataU(end).type = 'float';
    % strain xy
    scalardataU(end+1).name = 'exy_us';
    scalardataU(end).data = stressU(:,3);
    scalardataU(end).type = 'float';
    % stress xx
    scalardataU(end+1).name = 'Sxx_us';
    scalardataU(end).data = stressU(:,1);
    scalardataU(end).type = 'float';
    % stress yy
    scalardataU(end+1).name = 'Syy_us';
    scalardataU(end).data = stressU(:,2);
    scalardataU(end).type = 'float';
    % stress xy
    scalardataU(end+1).name = 'Sxy_us';
    scalardataU(end).data = stressU(:,3);
    scalardataU(end).type = 'float';
   
    % displacement in x (scalar)
    scalardataU(end+1).name = 'disp_usx';
    scalardataU(end).data = Solution.u(xdofs_u);
    scalardataU(end).type = 'float';
    % displacement in y (scalar)
    scalardataU(end+1).name = 'disp_usy';
    scalardataU(end).data = Solution.u(ydofs_u);
    scalardataU(end).type = 'float';
    % velocity in x (scalar)
    scalardataU(end+1).name = 'vel_usx';
    scalardataU(end).data = Solution.udot(xdofs_u);
    scalardataU(end).type = 'float';
    % velocity in y (scalar)
    scalardataU(end+1).name = 'vel_usy';
    scalardataU(end).data = Solution.udot(ydofs_u);
    scalardataU(end).type = 'float';
    
    %% Fluid displacement - 2D case
    if contains(Control.PMmodel, 'UPU')
        % fixed DOFs
        scalardataUf(1).name = 'fixed_uf';
        scalardataUf(1).data = fixedUf(MeshU.DOF);
        scalardataUf(1).type = 'int';
        % displacement
        vectordataUf(1).name = 'disp_uf';
        vectordataUf(1).data = [Solution.uf(xdofs_u) Solution.uf(ydofs_u) zeros(length(xdofs_u),1)];
        vectordataUf(1).type = 'float';
        % velocity
        vectordataUf(end+1).name = 'vel_uf';
        vectordataUf(end).data = [Solution.ufdot(xdofs_u) Solution.ufdot(ydofs_u) zeros(length(xdofs_u),1)];
        vectordataUf(end).type = 'float';
        
        if contains(Control.PMmodel, 'Dyn')
            % acceleration
            vectordataUf(end+1).name = 'acc_uf';
            vectordataUf(end).data = [Solution.uf2dot(xdofs_u) Solution.uf2dot(ydofs_u) zeros(length(xdofs_u),1)];
            vectordataUf(end).type = 'float';
        end
        
        % displacement in x (scalar)
        scalardataUf(end+1).name = 'disp_ufx';
        scalardataUf(end).data = Solution.uf(xdofs_u);
        scalardataUf(end).type = 'float';
        % displacement in y (scalar)
        scalardataUf(end+1).name = 'disp_ufy';
        scalardataUf(end).data = Solution.uf(ydofs_u);
        scalardataUf(end).type = 'float';
        % velocity in x (scalar)
        scalardataUf(end+1).name = 'vel_ufx';
        scalardataUf(end).data = Solution.ufdot(xdofs_u);
        scalardataUf(end).type = 'float';
        % velocity in y (scalar)
        scalardataUf(end+1).name = 'vel_ufy';
        scalardataUf(end).data = Solution.ufdot(ydofs_u);
        scalardataUf(end).type = 'float';        
    end 
    
    %% Fluid velocity - 2D case
    if contains(Control.PMmodel, 'UPV')
        % fixed DOFs
        scalardataUfdot(1).name = 'fixed_ufdot';
        scalardataUfdot(1).data = fixedUfdot(MeshU.DOF);
        scalardataUfdot(1).type = 'int';
        % velocity
        vectordataUfdot(1).name = 'vel_uf';
        vectordataUfdot(1).data = [Solution.ufdot(xdofs_u) Solution.ufdot(ydofs_u) zeros(length(xdofs_u),1)];
        vectordataUfdot(1).type = 'float';
        
        if contains(Control.PMmodel, 'Dyn')
            % acceleration
            vectordataUfdot(end+1).name = 'acc_uf';
            vectordataUfdot(end).data = [Solution.uf2dot(xdofs_u) Solution.uf2dot(ydofs_u) zeros(length(xdofs_u),1)];
            vectordataUfdot(end).type = 'float';
        end
        
        % velocity in x (scalar)
        scalardataUfdot(end+1).name = 'vel_ufx';
        scalardataUfdot(end).data = Solution.ufdot(xdofs_u);
        scalardataUfdot(end).type = 'float';
        % velocity in y (scalar)
        scalardataUfdot(end+1).name = 'vel_ufy';
        scalardataUfdot(end).data = Solution.ufdot(ydofs_u);
        scalardataUfdot(end).type = 'float';        
    end 
    
    %% Fluid velocity (relative) - 2D case
    if contains(Control.PMmodel, 'UPW')
        % fixed DOFs
        scalardataW(1).name = 'fixed_w';
        scalardataW(1).data = fixedW(MeshU.DOF);
        scalardataW(1).type = 'int';
        % velocity
        vectordataW(1).name = 'vel_w';
        vectordataW(1).data = [Solution.w(xdofs_u) Solution.w(ydofs_u) zeros(length(xdofs_u),1)];
        vectordataW(1).type = 'float';
        
        if contains(Control.PMmodel, 'Dyn')
            % acceleration
            vectordataW(end+1).name = 'acc_w';
            vectordataW(end).data = [Solution.wdot(xdofs_u) Solution.wdot(ydofs_u) zeros(length(xdofs_u),1)];
            vectordataW(end).type = 'float';
        end
        
        % velocity in x (scalar)
        scalardataW(end+1).name = 'vel_wx';
        scalardataW(end).data = Solution.w(xdofs_u);
        scalardataW(end).type = 'float';
        % velocity in y (scalar)
        scalardataW(end+1).name = 'vel_wy';
        scalardataW(end).data = Solution.w(ydofs_u);
        scalardataW(end).type = 'float';        
    end 
    
end

%% Storing data for fluid media
% pressure
scalardataP(1).name = 'pressure';
scalardataP(1).data = Solution.p;
scalardataP(1).type = 'float';
% fixed DOFs
scalardataP(end+1).name = 'fixedP';
scalardataP(end).data = fixedP(MeshP.DOF);
scalardataP(end).type = 'int';

if nsd == 1
    % fluid flux in 1D
    scalardataP(end+1).name = 'q';
    scalardataP(end).data = q;
    scalardataP(end).type = 'float';
    vectordataP = [];
elseif nsd == 2
    % fluid flux in 2D
    vectordataP(1).name = 'q';
    vectordataP(1).data = q;
    vectordataP(1).type = 'float';
end

%% Storing data for porosity field
% porosity
if contains(Control.PMmodel, 'UPN')
    scalardataN(1).name = 'eta';
    scalardataN(1).data = Solution.n;
    scalardataN(1).type = 'float';
elseif contains(Control.PMmodel, 'UPU') || contains(Control.PMmodel, 'UPV') || contains(Control.PMmodel, 'UPW')
    scalardataN(1).name = 'eta';
    scalardataN(1).data = eta;
    scalardataN(1).type = 'float';
    
    scalardataN(end+1).name = 'eta_dot';
    scalardataN(end).data = etadot;
    scalardataN(end).type = 'float';
else
    scalardataN = [];
end

%% Material type
scalardataU(end+1).name = 'mat_type';
scalardataU(end).data = MeshU.MatNodes;
scalardataU(end).type = 'int';

%% Write to VTK
% solid
description = config_name; % config file name
nameU = 'Solution_us.vtk.';
filenameU = fullfile(vtk_dir, [nameU, num2str(step)]);
WriteMesh2VTK(filenameU, description, MeshU, scalardataU, vectordataU);

% fluid
description = config_name; % config file name
nameP = 'Solution_p.vtk.';
filenameP = fullfile(vtk_dir, [nameP, num2str(step)]);
WriteMesh2VTK(filenameP, description, MeshP, scalardataP, vectordataP);

% porosity
if contains(Control.PMmodel, 'UPN') 
    description = config_name; % config file name
    nameN = 'Solution_n.vtk.';
    filenameN = fullfile(vtk_dir, [nameN, num2str(step)]);
    WriteMesh2VTK(filenameN, description, MeshN, scalardataN, []);
elseif contains(Control.PMmodel, 'UPU') || contains(Control.PMmodel, 'UPV') || contains(Control.PMmodel, 'UPW')
    description = config_name; % config file name
    nameN = 'Solution_n.vtk.';
    filenameN = fullfile(vtk_dir, [nameN, num2str(step)]);
    WriteMesh2VTK(filenameN, description, MeshU, scalardataN, []);
end

% fluid (u-p-u)
if contains(Control.PMmodel, 'UPU')
    description = config_name; % config file name
    nameUf = 'Solution_uf.vtk.';
    filenameUf = fullfile(vtk_dir, [nameUf, num2str(step)]);
    WriteMesh2VTK(filenameUf, description, MeshU, scalardataUf, vectordataUf);
end

% fluid (u-p-v)
if contains(Control.PMmodel, 'UPV')
    description = config_name; % config file name
    nameUfdot = 'Solution_ufdot.vtk.';
    filenameUfdot = fullfile(vtk_dir, [nameUfdot, num2str(step)]);
    WriteMesh2VTK(filenameUfdot, description, MeshU, scalardataUfdot, vectordataUfdot);
end

% fluid (u-p-w)
if contains(Control.PMmodel, 'UPW')
    description = config_name; % config file name
    nameW = 'Solution_w.vtk.';
    filenameW = fullfile(vtk_dir, [nameW, num2str(step)]);
    WriteMesh2VTK(filenameW, description, MeshU, scalardataW, vectordataW);
end

end