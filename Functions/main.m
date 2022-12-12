% Porous Media Simulation
% ------------------------------------------------------------------------
% Created by Bruna Campos
% bccampos@uwaterloo.ca
% Department of Civil Engineering, University of Waterloo
% January 2022
% ------------------------------------------------------------------------
% Reference: https://github.com/GCMLab
% ------------------------------------------------------------------------

%% Delete past files
if plot2vtk || plot2csv_on
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
[Material, MeshU, MeshP, BC, Control] = feval(config_name, ConfigDir, progress_on);

%% Quadrature points
Quad = GlobalQuad(MeshU, Control);

%% Assemble system matrices
disp([num2str(toc),': Assembling System Matrices...']);

% choice of quasi-steady or transient case
if Control.steady
    % quasi-steady case
    [Kuu, Kup, Kpp, S] = ComputeSystemMatrices_Steady(Material, MeshU, MeshP, Control, Quad);
else
    % transient case (mass matrices computed)
    [Kuu, Kup, Kpp, M, Mhat, S] = ComputeSystemMatrices_Dynamic(Material, MeshU, MeshP, Control, Quad);
end

% load vectors
[fu,fp] = ComputeSystemLoads(BC, MeshU, MeshP);

%% Initialize iteration variables
nsd = MeshU.nsd; % number of spatial directions

Iteration.u_old = zeros(MeshU.nDOF, 1);
Iteration.p_old = zeros(MeshP.nDOF, 1);

if ~Control.steady
    Iteration.udot_old = zeros(MeshU.nDOF, 1);
    Iteration.u2dot_old = zeros(MeshU.nDOF, 1);
    Iteration.pdot_old = zeros(MeshP.nDOF, 1);
end

% arrays for plot in time
Plot.time = (0:Control.dt:Control.tend);
Plot.p = zeros(length(Plot.time), 1); % fluid pressure
Plot.pan = zeros(length(Plot.time), 1); % analytic fluid pressure (quasi-steady case)
Plot.u = zeros(length(Plot.time), 1); % solid displacement
Plot.udot = zeros(length(Plot.time), 1); % solid velocity

%% Solution parameters
Control.t = 0;  % initial simulation time
Control.step = 1; % step counter

%% Solve system
while Control.t < Control.tend
    fprintf('\n Step %d \n', Control.step);

    % linear solver
    if Control.steady
        if MeshU.nsd == 1
            % analytical solution for 1D case
            [~, p_an, u_an] = getAnalyticResult (Material, MeshU, MeshP, BC, Control);
            Plot.p_an = p_an;
            Plot.u_an = u_an;
            Plot.pan(Control.step+1) = p_an(Control.plotp, 1);
        end
        % quasi-steady case
        [u,udot, p, pdot,~] = SolverSteady(Kuu, Kup, Kpp, S, fu, fp, BC, Control, Iteration);
    else
        % transient case
        [u, udot, u2dot, p, pdot, ~] = SolverDynamic(Kuu, Kup, Kpp, M, Mhat, S, fu, fp, BC, Control, Iteration);
    end

    % store variables over time
    if Control.step < length(Plot.time)
        % plot pressure vs time
        Plot.p(Control.step+1) = p(Control.plotp, 1);
        % plot displacement vs time
        Plot.u(Control.step+1) = u(Control.plotu, 1);
        % plot velocity vs time
        Plot.udot(Control.step+1) = udot(Control.plotu, 1);
    end

    % post processing: export VTK file
    if plot2vtk
        PostProcessing(u, udot, u2dot, p, Material, MeshU, MeshP, Control, Quad, BC, config_name, vtk_dir);
    end

    % update variables
    Iteration.u_old = u;
    Iteration.p_old = p;
    if ~Control.steady
        Iteration.udot_old = udot;
        Iteration.u2dot_old = u2dot;
        Iteration.pdot_old = pdot;
    end

    % update time and step
    Control.t = Control.t + Control.dt;
    Control.step = Control.step + 1;

end

%% Post processing

% plot figures
PlotGraphs(MeshU, MeshP, Control, Plot, u, p, saveGraphs_on);

% export CSV file
if plot2csv_on
    plot2csv(config_name, vtk_dir, Plot);
end