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
[Material, MeshU, MeshP, MeshN, BC, Control] = feval(config_name, ConfigDir, progress_on);

%% Quadrature points
Quad = GlobalQuad(MeshU, Control);

%% Solution parameters
Control.t = 0;  % initial simulation time
Control.step = 0; % step counter

%% Assemble system matrices
disp([num2str(toc),': Assembling System Matrices...']);

if Control.steady %% quasi-steady/transient case
    if Control.Biotmodel
        [Kuu, Kup, Kpp, S] = ComputeSystemMatrices_Transient(Material, MeshU, MeshP, Quad);
    else
        [Kuu, Kup, Kpp, Kpu, S, Kpn, Knn, Knu, Knp, Aun, Bun] = ComputeSystemMatrices_Transient_v2(Material, MeshU, MeshP, MeshN, Quad);
    end
else %% dynamic case
    [Kuu, Kup, Kpp, M, Mhat, S] = ComputeSystemMatrices_Dynamic(Material, MeshU, MeshP, Quad);
end

%% Assemble system load vectors
[fu,fp,fn] = ComputeSystemLoads(BC, MeshU, MeshP, MeshN, Control, Quad);

%% Initialize iteration variables
nsd = MeshU.nsd; % number of spatial directions

Iteration.u_old = zeros(MeshU.nDOF, 1); % displacement variable storage
Iteration.udot_old = zeros(MeshU.nDOF, 1); % solid velocity
Iteration.p_old = zeros(MeshP.nDOF, 1); % pressure variable storage
Iteration.u2dot_old = zeros(MeshU.nDOF, 1); % solid acceleration
Iteration.pdot_old = zeros(MeshP.nDOF, 1); % pressure gradient

% arrays for plot in time
Plot.time = (0:Control.dt:Control.tend);
Plot.p = zeros(length(Plot.time), 1); % fluid pressure
Plot.pan = zeros(length(Plot.time), 1); % analytic fluid pressure (quasi-steady case)
Plot.u = zeros(length(Plot.time), 1); % solid displacement
Plot.udot = zeros(length(Plot.time), 1); % solid velocity

%% Initial conditions
% displacement
if ~isempty(BC.initU)
    Iteration.u_old = BC.initU;
    Plot.u = BC.initU(Control.plotu,1);
end

% pressure
if ~isempty(BC.initP)
    Iteration.p_old = BC.initP;
    Plot.p(1,1) = BC.initP(Control.plotp,1);
end

% porosity
if ~Control.Biotmodel
    % initial porosity condition
    Iteration.n_old = zeros(MeshN.nDOF, 1);
    Iteration.n_old(:) = Material.n;
    % plot porosity over time
    Plot.n = zeros(length(Plot.time), 1);
    Plot.n(1,1) = Material.n;
end

%% Initial condition file
if plot2vtk
    Solution.u = Iteration.u_old;
    Solution.udot = Iteration.udot_old;
    Solution.u2dot = Iteration.u2dot_old;
    Solution.p = Iteration.p_old;
    if ~Control.Biotmodel
        Solution.n = Iteration.n_old;
    else
        Solution.n = [];
    end
    % save vtk file
    PostProcessing(Solution, Material, MeshU, MeshP, MeshN, Control, BC, config_name, vtk_dir);
    % update time step
    Control.step = 1;
end
    
%% Solve system
while Control.t < Control.tend
    fprintf('\n Step %d \n', Control.step);
    
    %%%% change for this load computation if time dependent
%     [fu,fp,fn] = ComputeSystemLoads(BC, MeshU, MeshP, MeshN, Control, Quad);
    
    % linear solver
    if Control.steady
%         if MeshU.nsd == 1 && Material.Qinv == 0
%             % analytical solution for 1D case and 1/Q = 0 (null S matrix)
%             % (incompressible solid and fluid)
%             [~, p_an, u_an] = getAnalyticResult (Material, MeshU, MeshP, BC, Control);
%             Plot.p_an = p_an;
%             Plot.u_an = u_an;
%             % store variables over time
%             if Control.step < length(Plot.time)
%                 Plot.pan(Control.step+1) = p_an(Control.plotp, 1);
%             end
%         end
        % quasi-steady/transient case
        if Control.Biotmodel
            [Solution] = SolverTransient(Kuu, Kup, Kpp, S, fu, fp, BC, Control, Iteration);
        else
            [Solution] = SolverTransient_v2(Kuu, Kup, Kpp, Kpu, S, Kpn, Knn, Knu, Knp, Aun, Bun, fu, fp, fn, BC, Control, Iteration);
        end
    else
        % dynamic case
        [Solution] = SolverDynamic_v4(Kuu, Kup, Kpp, M, Mhat, S, fu, fp, BC, Control, Iteration);
    end
    
    % update external forces vectors
    fu(BC.fixed_u) = Solution.fE;
    fp(BC.fixed_p) = Solution.qE;
    
    % store variables over time
    if Control.step < length(Plot.time)
        % plot pressure vs time
        Plot.p(Control.step+1) = Solution.p(Control.plotp, 1);
        % plot displacement vs time
        Plot.u(Control.step+1) = Solution.u(Control.plotu, 1);
        % plot velocity vs time
        Plot.udot(Control.step+1) = Solution.udot(Control.plotu, 1);
        % plot porosity vs time
        if ~Control.Biotmodel
            Plot.n(Control.step+1) = Solution.n(Control.plotp, 1);
        end
    end
    
    % post processing: compute stress/flux, export VTK file
    if plot2vtk
        PostProcessing(Solution, Material, MeshU, MeshP, MeshN, Control, BC, config_name, vtk_dir);
    end
    
    % update variables
    Iteration.u_old = Solution.u; % solid displacement
    Iteration.udot_old = Solution.udot; % solid velocity
    Iteration.p_old = Solution.p; % fluid pressure
    if ~Control.Biotmodel
        Iteration.n_old = Solution.n; % medium porosity
    end
    if ~Control.steady
        Iteration.u2dot_old = Solution.u2dot; % solid acceleration
        Iteration.pdot_old = Solution.pdot; % fluid pressure gradient
    end
    
    % update time and step
    Control.t = Control.t + Control.dt;
    Control.step = Control.step + 1;
    
end

%% Post processing

% plot figures
PlotGraphs(Solution, MeshU, MeshP, MeshN, Control, Plot, saveGraphs_on);

% export CSV file
if plot2csv_on
    Plot2csv(config_name, vtk_dir, Plot);
end
