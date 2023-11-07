% Porous Media Simulation
% ------------------------------------------------------------------------
% Created by Bruna Campos
% bccampos@uwaterloo.ca
% Department of Civil Engineering, University of Waterloo
% January 2022
% ------------------------------------------------------------------------
% Reference: https://github.com/GCMLab
% ------------------------------------------------------------------------
% Version 6 (07/26/22): update porosity field for time and spatial
% variation

%% Clear variables and initialize code
clearvars
clear, clc
close all
format compact
tic;

% use current directory
curDir = pwd;

%% Input data
% Config files folder
DirFolder = 'Config Files';
% Config file to run
% -------------------- Consolidation 1D - quasi-steady
% File = 'Column1D_Steady_Boone';
% File = 'Column1D_Steady_Ferronato';
% File = 'Column1D_Steady_Korsawe';
% File = 'Column1D_Steady_Tasiop';

% -------------------- Consolidation 1D - dynamic
% File = 'Column1D_Dynamic_Boone';
% File = 'Column1D_Dynamic_Diebels';
% File = 'Column1D_Dynamic_Gajo';
% File = 'Column1D_Dynamic_Komijani';
% File = 'Column1D_Dynamic_Sandstone';
% File = 'Column1D_Dynamic_Tasiop';
% File = 'Column1D_Dynamic_Zienck';

% -------------------- Consolidation 2D - quasi-steady
% File = 'Column2D_Steady_Ferronato';
% File = 'Column2D_Steady_Zheng';

% -------------------- Consolidation 2D - dynamic
% File = 'Column2D_Dynamic_Komijani';

% -------------------- Footing
% File = 'Footing2D_Diebels';
% File = 'Footing2D_Korsawe';

% -------------------- Injection well
% File = 'Injection2D_Plate';
% File = 'InjectionWellPoint_Dynamic_Komijani';
% File = 'InjectionWells2D_15x15';
% File = 'Plate_Dynamic_Komijani';
% File = 'PlatePointInjection_Dynamic';
% File = 'PlatePointInjection_Dynamic_Komijani';

% -------------------- Tests Convergence
% File = 'ManufacturedSolution1D_Biot';
% File = 'ManufacturedSolution1D_Pste';
% File = 'ManufacturedSolution1D_UstePste';
% File = 'ManufacturedSolution1D_UstePtra';
% File = 'ManufacturedSolution1D_UtraPtra';
% File = 'ManufacturedSolution1Dupu_Biot';
% File = 'ManufacturedSolution1Dupu_Spanos';

% -------------------- Tests Elasticity
% File = 'Beam_Dynamic';
% File = 'Beam_Dynamic_v2';
% File = 'Plate2D_KirschTest';
% File = 'PlateWithHole_Elasticity';
% File = 'VelocityImpact1D_Ham';
% File = 'VelocityImpact1D_Idesman';

% -------------------- Tests Heat Transfer
% File = 'HeatConduction1D_Dynamic';
% File = 'HeatConduction2D';
% File = 'HeatConduction2D_PointFluxBC';
% File = 'HeatConduction2D_TempBC';
% File = 'PlateDiffusionDynamic';
% File = 'PlateDiffusionDynamic_Point';
% File = 'PlateDiffusionSteady';
% File = 'PlateWithHole_Diffusion';
% File = 'PlateWithHole_Diffusion_Transient';
% File = 'PlateWithHole_HeatTransfer';

% -------------------- Velocity Impact
% File = 'VelocityImpact2D_Dynamic_Komijani';
% File = 'VelocityImpact2D_Dynamic_Plate';

% -------------------- Wave Propagation
% File = 'Column1D_Dynamic_Pulse';
% File = 'Column1D_Dynamic_Pulse5000m';
% File = 'Plate2D_Dynamic_Pulse';
% File = 'WaveProp_Dynamic_Komijani';
% File = 'WaveProp_Dynamic_Quiroga';
% File = 'WaveProp_Dynamic_Tian';
File = 'WaveProp_Dynamic_InjectionPress10m';
% File = 'WaveProp_Dynamic_InjectionPress100m';
% File = 'WaveProp_Dynamic_InjectionFlux';

% ------------------------------------------------------------------------

% Directory for VTK file
% VTKFolder ='C:\Users\bu_ca\Downloads\PoroelasticityFEM\Results';
VTKFolder ='C:\Users\bccampos\Downloads\PoroelasticityFEM\Results';

% output VTK file
plot2vtk = 1;
% export Matlab images
saveGraphs_on = 0;
% export Matlab files
saveMatData_on = 0;
% output progress messages
progress_on = 1;
% export video
saveVideo_on = 0;

%% Directories
VTKFolder = fullfile(VTKFolder, DirFolder);
FuncDir = fullfile(curDir, 'Functions');
ConfigDir = fullfile(curDir, DirFolder);

%% Add validation folder and function folder to the search path
% genpath adds all subfolders as well
addpath(genpath(FuncDir));
addpath(genpath(ConfigDir));

ConfigFile = {File};
VTKDir = fullfile(VTKFolder,ConfigFile);

%% Run config file
% clear variables from previous simulation
clearvars -except VTKDir ConfigFile...
    curDir FuncDir  ConfigDir ...
    file codeSubmitTime ...
    exit_when_done print_log ...
    plot2vtk progress_on ...
    saveGraphs_on saveMatData_on ...
    saveVideo_on

clearvars -global

% file name
config_name_full = ConfigFile{1};
[~,config_name] = fileparts(config_name_full);

% create VTK folder for config file, if not existent
vtk_dir = VTKDir{1};
if ~isfolder(vtk_dir)
    mkdir(vtk_dir)
end

% run and time the simulation
start_time = toc;
run('Functions/Main/main');
end_time = toc;

disp(['run time: ' num2str(end_time - start_time)])
% close all