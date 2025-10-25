% This file is part of PoroelasticityFEM.
%
% PoroelasticityFEM is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% PoroelasticityFEM is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% Copyright (C) 2022-2024 Bruna Campos

% ------------------------------------------------------------------------
% Porous Media Simulation - Tests
% ------------------------------------------------------------------------
% Created by Bruna Campos
% bccampos@uwaterloo.ca
% Department of Civil Engineering, University of Waterloo
% January 2022
% ------------------------------------------------------------------------
% Reference: https://github.com/GCMLab
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
VTKFolder = pwd;

% output VTK file
plot2vtk = 0;
% export Matlab images
saveGraphs_on = 0;
% export Matlab files
saveMatData_on = 0;
% output progress messages
progress_on = 1;
% plot graphs
plotGraphs_on = 0;
% save vide file
saveVideo_on = 0;

%% Directories
VTKFolder = fullfile(VTKFolder, DirFolder);
FuncDir = fullfile(curDir, 'Functions');
ConfigDir = fullfile(curDir, DirFolder);

%% Add validation folder and function folder to the search path
% genpath adds all subfolders as well
addpath(genpath(FuncDir));
addpath(genpath(ConfigDir));

% number of tests
ntests = 12;
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
run('Test Files/RunT6_PatchTestF')

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

%% Test 11: Manufactured Solution - Transient case convergence
%   Pass condition: L2-norm converges at a rate of at least h^3
%                    e-norm converges at a rate of at least h^2
run('Test Files/RunT11_ManSolTransient')

%% Test 12: Poroelasticity - Consolidation problem
%   Pass condition: Displacement and pressure are the same compared to
%   analytical solution
run('Test Files/RunT12_PMConsolidation')

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