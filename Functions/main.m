% Porous Media Simulation
% Bruna Campos
% January 2022
%
% Main Function
% Adapted from: https://github.com/GCMLab (Acknowledgements: Matin Parchei
% Esfahani)

%% Delete past files
if plot2vtk || plot2csv_on
    if progress_on
        fprintf('%.2f: Deleting old vtk files...\n', toc);
    end
    sol = fullfile(vtk_dir, '*');
    delete(sol)
end

%% Input problem
if progress_on
    fprintf('%.2f: Reading config file...\n', toc);
end
[Material, Mesh, BC, Control] = feval(config_name);

%% Quadrature points
Quad = GlobalQuad(Control);

%% Assemble system matrices
disp([num2str(toc),': Assembling System Matrices...']);

% choice of quasi-steady or transient case
if Control.steady
    % quasi-steady case
    [Kuu, Kup, Kpp, S] = ComputeSystemMatrices_Steady(Material, Mesh, Control, Quad);
else
    % transient case (mass matrices computed)
    [Kuu, Kup, Kpp, M, Mhat, S] = ComputeSystemMatrices_Dynamic(Material, Mesh, Control, Quad);
end

%% Initialize iteration variables
Iteration.u_old = zeros(Mesh.ndof_u, 1);
Iteration.p_old = zeros(Mesh.ndof_p,1);

if ~Control.steady
    Iteration.udot_old = zeros(Mesh.ndof_u,1);
    Iteration.u2dot_old = zeros(Mesh.ndof_u,1);
    Iteration.pdot_old = zeros(Mesh.ndof_p,1);
end

% arrays for plot in time
Plot.time = (0:Control.dt:Control.tend);
Plot.p = zeros(1, length(Plot.time)); % fluid pressure
Plot.pan = zeros(1, length(Plot.time)); % analytic fluid pressure (quasi-steady case)
Plot.u = zeros(1, length(Plot.time)); % solid displacement
Plot.udot = zeros(1, length(Plot.time)); % solid velocity

%% Solution parameters
Control.t = 0;  % initial simulation time
Control.step = 1; % step counter

%% Solve system
while Control.t < Control.tend

    fprintf('\n Step %d \n', Control.step);

    % linear solver
    if Control.steady
        % analytical solution
        [~, p_an, u_an] = getAnalyticResult (Material, Mesh, BC, Control);
        Plot.p_an = p_an;
        Plot.u_an = u_an;
        Plot.pan(Control.step+1) = p_an(1, Control.plotp);
        
        % quasi-steady case
        [u,udot, p, pdot,~] = SolverSteady(Kuu, Kup, Kpp, S, Mesh, BC, Control, Iteration);
    else
        % transient case
        [u, udot, u2dot, p, pdot, ~] = SolverDynamic(Kuu, Kup, Kpp, M, Mhat, S, Mesh, BC, Control, Iteration);

        % TEST
%         [u, udot, u2dot, p, pdot, ~] = SolverFrequency(Kuu, Kup, Kpp, M, Mhat, S, Mesh, BC, Control, Iteration);
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
postprocessor(Mesh, Control, Plot, u, p);

% export CSV file
if plot2csv_on
    plot2csv(config_name, vtk_dir, Plot);
end