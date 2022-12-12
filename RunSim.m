% Porous Media Simulation
% ------------------------------------------------------------------------
% Created by Bruna Campos
% bccampos@uwaterloo.ca
% Department of Civil Engineering, University of Waterloo
% January 2022
% ------------------------------------------------------------------------
% Reference: https://github.com/GCMLab
% ------------------------------------------------------------------------
% Version 5 (July 2022): 1D Spanos quasi static model

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
% File = 'Column1D_Steady_Korsawe';
File = 'Column1D_Steady_Boone';
% File = 'Column1D_Dynamic_Komijani';
% File = 'Column1D_Dynamic_Zienck';
% File = 'Column2D_Dynamic_Komijani';

% Directory for VTK file
% VTKFolder ='C:\Users\bu_ca\OneDrive\Documents\Doutorado UWaterloo\Research\Matlab Codes\PorousMedia\Results\';
VTKFolder ='C:\Users\bccampos\OneDrive\Documents\Doutorado UWaterloo\Research\Matlab Codes\PorousMedia_v5\Results\';

% output VTK file
plot2vtk = 0;
% output CSV file
plot2csv_on = 0;
% export Matlab images
saveGraphs_on = 0;
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
    plot2vtk progress_on plot2csv_on ...
    saveGraphs_on

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
run('Functions/main');
end_time = toc;

disp(['run time: ' num2str(end_time - start_time)])
% close all