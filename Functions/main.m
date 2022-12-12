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
% [Material, MeshU, MeshP, BC, Control] = feval(config_name, ConfigDir, progress_on);
[Material, MeshU, MeshP, MeshN, BC, Control] = feval(config_name, ConfigDir, progress_on);

%% Quadrature points
Quad = GlobalQuad(MeshU, Control);

%% Assemble system matrices
disp([num2str(toc),': Assembling System Matrices...']);

if Control.steady %% quasi-steady case
    if Control.Biotmodel
        [Kuu, Kup, Kpp, S] = ComputeSystemMatrices_Steady(Material, MeshU, MeshP, Control, Quad);
    else
        [Kuu, Kup, Kpp, Kpu, S, Kpn, Knn, Knu, Knp, Kun] = ComputeSystemMatrices_Steady_v2(Material, MeshU, MeshP, MeshN, Control, Quad);
    end
else %% transient case
    [Kuu, Kup, Kpp, M, Mhat, S] = ComputeSystemMatrices_Dynamic(Material, MeshU, MeshP, Control, Quad);
end

%% Assemble system load vectors
[fu,fp] = ComputeSystemLoads(BC, MeshU, MeshP, Control, Quad);

%% Initialize iteration variables
nsd = MeshU.nsd; % number of spatial directions

Iteration.u_old = zeros(MeshU.nDOF, 1); % displacement variable storage
Iteration.p_old = zeros(MeshP.nDOF, 1); % pressure variable storage

if ~Control.steady
    Iteration.udot_old = zeros(MeshU.nDOF, 1); % solid velocity
    Iteration.u2dot_old = zeros(MeshU.nDOF, 1); % solid acceleration
    Iteration.pdot_old = zeros(MeshP.nDOF, 1); % pressure gradient
end

% arrays for plot in time
Plot.time = (0:Control.dt:Control.tend);
Plot.p = zeros(length(Plot.time), 1); % fluid pressure
Plot.pan = zeros(length(Plot.time), 1); % analytic fluid pressure (quasi-steady case)
Plot.u = zeros(length(Plot.time), 1); % solid displacement
Plot.udot = zeros(length(Plot.time), 1); % solid velocity

%% Porosity variables
if ~Control.Biotmodel
    % initial porosity condition
    Iteration.n_old = zeros(MeshN.nDOF, 1);
    Iteration.n_old(:,1) = Material.n; 
    % plot porosity over time
    Plot.n = zeros(length(Plot.time), 1);
    Plot.n(1,1) = Material.n; 
end

%% Solution parameters
Control.t = 0;  % initial simulation time
Control.step = 1; % step counter

%% Solve system
while Control.t < Control.tend
    fprintf('\n Step %d \n', Control.step);
    
    % linear solver
    if Control.steady
        if MeshU.nsd == 1
            % analytical solution for 1D case and 1/Q = 0 (null S matrix)
            % (incompressible solid and fluid)
            [~, p_an, u_an] = getAnalyticResult (Material, MeshU, MeshP, BC, Control);
            Plot.p_an = p_an;
            Plot.u_an = u_an;
            % store variables over time
            if Control.step < length(Plot.time)
                Plot.pan(Control.step+1) = p_an(Control.plotp, 1);
            end
        end
        % quasi-steady case
        if Control.Biotmodel
            [u,udot, p, pdot, fE,~] = SolverSteady(Kuu, Kup, Kpp, S, fu, fp, BC, Control, Iteration);
            fu(BC.fixed_u) = fE;
        else
            [u,udot, p, pdot, n, ndot, fE, ~] = SolverSteady_v2(Kuu, Kup, Kpp, S, Kpu, Kun, Kpn, Knn, Knu, Knp, fu, fp, BC, Control, Iteration);
            fu(BC.fixed_u) = fE;
        end
    else
        % transient case
        [u, udot, u2dot, p, pdot, ~, ~] = SolverDynamic(Kuu, Kup, Kpp, M, Mhat, S, fu, fp, BC, Control, Iteration);
    end
    
    % store variables over time
    if Control.step < length(Plot.time)
        % plot pressure vs time
        Plot.p(Control.step+1) = p(Control.plotp, 1);
        % plot displacement vs time
        Plot.u(Control.step+1) = u(Control.plotu, 1);
        % plot velocity vs time
        Plot.udot(Control.step+1) = udot(Control.plotu, 1);
        % plot porosity vs time
        if ~Control.Biotmodel
            Plot.n(Control.step+1) = n(Control.plotp, 1);
        end
    end
    
    % post processing: export VTK file
    if plot2vtk
        if Control.Biotmodel
            PostProcessing(u, udot, u2dot, p, Material, MeshU, MeshP, Control, Quad, BC, config_name, vtk_dir);
        else
            %%%% NEEDS TO BE ADAPTED FOR POROSITY
        end
    end
    
    % update variables
    Iteration.u_old = u;
    Iteration.p_old = p;
    if ~Control.Biotmodel
        Iteration.n_old = n;
    end
    
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
PlotGraphs(MeshU, MeshP, MeshN, Control, Plot, u, p, saveGraphs_on);

% export CSV file
if plot2csv_on
    plot2csv(config_name, vtk_dir, Plot);
end