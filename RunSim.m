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
% -------------------- 1D steady/transient column consolidation
% File = 'Column1D_Steady_Korsawe';
File = 'Column1D_Steady_Boone';
% File = 'Column1D_Steady_Quiroga';
% File = 'Column1D_Steady_Sandstone';
% File = 'Column1D_Steady_Zheng';
% File = 'Column1D_Steady_Ferronato';
% File = 'Column1D_Steady_Tasiop';

% -------------------- 1D dynamic column consolidation
% File = 'Column1D_Dynamic_Komijani';
% File = 'Column1D_Dynamic_Gajo';
% File = 'Column1D_Dynamic_Tasiop';
% File = 'Column1D_Dynamic_Zienck';
% File = 'Column1D_Dynamic_Diebels';

% -------------------- 1D dynamic velocity impact
% File = 'VelocityImpact1D_Ham';
% File = 'VelocityImpact1D_Idesman';
% File = 'VelocityImpact1D_Dynamic_Komijani';

% -------------------- 1D harmonic oscillation
% File = 'Harmonic1D';

% -------------------- 1D tests
% File = 'HeatConduction_1D_Transient';
% File = 'ManufacturedSolution1D_Biot';
% File = 'ManufacturedSolution1D_Pste';
% File = 'ManufacturedSolution1D_UstePste';
% File = 'ManufacturedSolution1D_UstePtra';
% File = 'ManufacturedSolution1D_UtraPtra';

% -------------------- 2D steady/transient column consolidation
% File = 'Column2D_Steady_Korsawe';
% File = 'Column2D_Steady_Boone';
% File = 'Column2D_Steady_Zheng';
% File = 'Column2D_Steady_Ferronato';

% -------------------- 2D dynamic column consolidation
% File = 'Column2D_Dynamic_Komijani';
% File = 'Plate_Dynamic_Komijani';

% -------------------- 2D injection wells
% File = 'InjectionWells2D_15x15';
% File = 'InjectionWellPoint_Dynamic_Komijani';
% File = 'PlatePointInjection_Dynamic_Komijani';
% File = 'PlatePointInjection_Dynamic';

% -------------------- 2D footing
% File = 'Footing2D_Korsawe';
% File = 'Footing2D_Diebels';

% -------------------- 2D dynamic velocity impact
% File = 'VelocityImpact2D_Dynamic_Komijani';
% File = 'VelocityImpact2D_Dynamic_Plate';

% -------------------- 2D wave propagation
% File = 'WaveProp_Dynamic_Quiroga';
% File = 'WaveProp_Dynamic_Komijani';

% -------------------- 2D tests: elasticity
% File = 'PlateWithHole_Elasticity';
% File = 'PlateWithHole_HeatTransfer';
% File = 'Beam_Dynamic';
% File = 'Beam_Dynamic_v2';
% File = 'PlateWithHole_Elasticity_Stress';
% File = 'Plate2D_KirschTest';

% -------------------- 2D tests: diffusion
% File = 'PlateWithHole_Diffusion';
% File = 'PlateWithHole_Diffusion_Transient';
% File = 'PlateDiffusionSteady';
% File = 'PlateDiffusionDynamic';
% File = 'PlateDiffusionDynamic_Point';
% File = 'HeatConduction1D_Dynamic';
% File = 'HeatConduction2D';
% File = 'HeatConduction2D_TempBC';
% File = 'HeatConduction2D_PointFluxBC';

% Directory for VTK file
VTKFolder ='C:\Users\bu_ca\Downloads\PoroelasticityFEM\Results';
% VTKFolder ='C:\Users\bccampos\Downloads\PoroelasticityFEM\Results';

% output VTK file
plot2vtk = 1;
% export Matlab images
saveGraphs_on = 0;
% export Matlab files
saveMatData_on = 1;
% output progress messages
progress_on = 1;

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
    saveGraphs_on saveMatData_on

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