% de la Cruz and Spanos poroelasticity model
% Dynamic case
% December 2022
% ------------------------------------------------------------------------
% May 2023: version with u-p system, substituting porosity equation
% ------------------------------------------------------------------------

%% Model name and type
disp([num2str(toc),': Model: Spanos u-p dynamic case']);

%% Assemble system matrices
disp([num2str(toc),': Assembling System Matrices...']);
[Muu, Mpu, Kuu, Kup, Kpp, Kpu, S] = ComputeMatricesDyn3_Spanos_UP(Material, MeshU, MeshP, QuadU);

%% Initial condition
% print initial condition step
fprintf('\n Step 1, t = 0 \n');

% initialize variables (account for initial condition)
[Iteration, Plot] = initVariables(MeshU, MeshP, MeshN, Material, Control, BC);

if plot2vtk
    % time domain solution
    Solution.u = Iteration.u_old;
    Solution.udot = Iteration.udot_old;
    Solution.u2dot = Iteration.u2dot_old;
    Solution.p = Iteration.p_old;

    % save vtk file
    PostProcessing(Solution, Material, MeshU, MeshP, MeshN, Control, BC, config_name, vtk_dir);
end

% initialize video file
if saveVideo_on
    myVideo = VideoWriter('myVideoFile'); %open video file
    myVideo.FrameRate = 20;
    open(myVideo)
end

%% Solve system
% start at t=dt (t=0 covered in initial condition)
for t = 2:length(Plot.time)
    % current time
    Control.t = Plot.time(t);
    % update time step
    Control.step = Control.step + 1;
    % print current time and step
    fprintf('\n Step %d, t = %d \n', Control.step, Control.t);
    
    % adapt time step size (dtc = dt current)
    if isfield(Control, 'dtmin') 
        if Control.t < Control.tlim
            Control.dtc = Control.dtmin;
        else
            Control.dtc = Control.dt;
        end
    else
        Control.dtc = Control.dt;
    end

    % system load vectors
    [fu,fp,fn] = ComputeLoads(BC, MeshU, MeshP, MeshN, Control, QuadU, QuadP);

    % linear solver
    [Solution] = SolverDyn_UP(Kuu, Kup, Kpp, Muu, Mpu, S, fu, fp, BC, Control, Iteration);

    % plot solution over time
    PlotGraphsUP_OverTime(MeshU, MeshP, Control, Solution);
    pause(0.0001);

    % update external forces vectors
    fu(BC.fixed_u) = Solution.fE;
    fp(BC.fixed_p) = Solution.qE;

    % post processing: compute stress/flux, export VTK file
    if plot2vtk
        PostProcessing(Solution, Material, MeshU, MeshP, MeshN, Control, BC, config_name, vtk_dir);
    end

    % store pressure over time
    Plot.p_time(Control.step,:) = Solution.p(Control.plotp, 1);
    % store displacement over time
    Plot.u_time(Control.step,:) = Solution.u(Control.plotu, 1);
    % store velocity over time
    Plot.udot_time(Control.step,:) = Solution.udot(Control.plotu, 1);
    % store solid acceleration over time
    Plot.u2dot_time(Control.step,:) = Solution.u2dot(Control.plotu, 1);

    % store variables over time at fixed row
    Plot.u_synthetic(Control.step,:) = Solution.u(Control.ploturow);
    Plot.udot_synthetic(Control.step,:) = Solution.udot(Control.ploturow);
    Plot.p_synthetic(Control.step,:) = Solution.p(Control.plotprow);

    % store variables over space
    Plot.p_space = Solution.p(Control.plotp);
    Plot.u_space = Solution.u(Control.plotu);
    Plot.udot_space = Solution.udot(Control.plotu);
    Plot.u2dot_space = Solution.u2dot(Control.plotu);

    % update variables
    Iteration.u_old = Solution.u; % solid displacement
    Iteration.udot_old = Solution.udot; % solid velocity
    Iteration.p_old = Solution.p; % fluid pressure
    Iteration.u2dot_old = Solution.u2dot; % solid acceleration
    Iteration.pdot_old = Solution.pdot; % fluid pressure gradient
    Iteration.fu_old = fu; % load vector
    Iteration.fp_old = fp; % flux vector
end

% store variables at fixed row
Plot.urow = Solution.u(Control.ploturow);
Plot.udotrow = Solution.udot(Control.ploturow);
Plot.prow = Solution.p(Control.plotprow);

% close video file
if saveVideo_on
    close(myVideo)
end
