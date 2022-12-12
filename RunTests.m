% Porous Media Simulation - Tests
% ------------------------------------------------------------------------
% Created by Bruna Campos
% bccampos@uwaterloo.ca
% Department of Civil Engineering, University of Waterloo
% January 2022
% ------------------------------------------------------------------------
% Reference: https://github.com/GCMLab (Acknowledgements: Bruce Gee)
% ------------------------------------------------------------------------

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
DirFolder = 'Test Files';

% Directory for VTK file
% VTKFolder ='C:\Users\bu_ca\OneDrive\Documents\Doutorado UWaterloo\Research\Poroelasticity codes\PorousMedia_v6\Results\';
VTKFolder ='C:\Users\bccampos\OneDrive\Documents\Doutorado UWaterloo\Research\Matlab Codes\PorousMedia_v6\Results\';

% output VTK file
plot2vtk = 1;
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

% number of tests
ntests = 10;
% initialize test summary
testpasssummary = zeros(ntests,1);
    
%% Run Tests
%% Test 1: patch test A
% Pass Condition: FEA solution displacements, stresses, and strains are exact
run('Test Files/RunT1_PatchTestA')

%% Test 2: patch test B
% Pass Condition: FEA solution displacements, stresses, and strains are exact
run('Test Files/RunT2_PatchTestB')

%% Test 3: patch test C
% Pass Condition: FEA solution displacements, stresses, and strains are exact
run('Test Files/RunT3_PatchTestC')

%% Test 4: patch test D
% Pass Condition: FEA solution displacements, stresses, and strains are exact
run('Test Files/RunT4_PatchTestD')

%% Test 5: patch test E
% Pass Condition: FEA solution displacements, stresses, and strains are exact
run('Test Files/RunT5_PatchTestE')

%% Test 6: patch test F
% Pass Condition: FEA solution displacements, stresses, and strains are exact
% run('Test Files/RunT6_PatchTestF')

%% Test 7: Manufactured Solution - Q4 element convergence
%   Pass condition: L2-norm converges at a rate of at least h^2
%                    e-norm converges at a rate of at least h
run('Test Files/RunT7_ManSolQ4')

%% Test 8: Manufactured Solution - Q9 element convergence
%   Pass condition: L2-norm converges at a rate of at least h^3
%                    e-norm converges at a rate of at least h^2
run('Test Files/RunT8_ManSolQ9')

%% Test 9: Manufactured Solution - L2 element convergence
%   Pass condition: L2-norm converges at a rate of at least h^2
%                    e-norm converges at a rate of at least h
run('Test Files/RunT9_ManSolL2')

%% Test 10: Manufactured Solution - L3 element convergence
%   Pass condition: L2-norm converges at a rate of at least h^3
%                    e-norm converges at a rate of at least h^2
run('Test Files/RunT10_ManSolL3')

%% Summarize test results
fprintf('\n\n%10s%10s', 'Test', 'Status')
fprintf('\n----------------------------------------------')
for test = 1:ntests
   if testpasssummary(test)
       fprintf('\n%10d%10s', test, 'PASS')
   else
      fprintf('\n%10d%10s', test, 'FAIL')
   end
end
fprintf('\n\n')