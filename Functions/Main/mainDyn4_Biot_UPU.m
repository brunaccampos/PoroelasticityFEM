% Biot poroelasticity model
% Dynamic case
% December 2022
% ------------------------------------------------------------------------

%% Model name and type
disp([num2str(toc),': Model: Biot u-p-U dynamic case']);

%% Assemble system matrices
disp([num2str(toc),': Assembling System Matrices...']);
[Kss, Ksp, Mss, Csf, Css, Kpf, Kps, Kpp, Kfp, Mff, Cff, Cfs, Msf, Mfs] = ComputeMatricesDyn4_Biot_UPU(Material, MeshU, MeshP, QuadU, QuadP);

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
    Solution.uf = Iteration.uf_old;
    Solution.ufdot = Iteration.ufdot_old;
    Solution.uf2dot = Iteration.uf2dot_old;
    
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
    [fu, fp, ff] = ComputeLoads_UPU(BC, MeshU, MeshP, Control, Material, QuadU, QuadP);
    
    % analytical solution for 1D case
    if Control.plotansol
        [p_an, u_an, uf_an] = feval(Control.ansol_type, Control, Material, MeshU, MeshP, BC);
       
        % store variables over length
        Plot.pan_space = p_an;
        Plot.uan_space = u_an;
        Plot.ufan_space = uf_an;
        
        % store variables over time
        Plot.uan_time(Control.step,:) = u_an(Control.plotu, 1);
        Plot.pan_time(Control.step,:) = p_an(Control.plotp, 1);
        Plot.ufan_time(Control.step,:) = uf_an(Control.plotu, 1);
    end
    
    % linear solver
    [Solution] = SolverDyn_UPU(Kss, Ksp, Mss, Csf, Css, Kpf, Kps, Kpp, Kfp, Mff, Cff, Cfs, Msf, Mfs, fu, fp, ff, BC, Control, Iteration);

    % plot solution over time
    PlotGraphsUPU_OverTime(MeshU, MeshP, Control, Solution);
    pause(0.0001);

    if saveVideo_on
        frame = getframe(gcf); %get frame
        writeVideo(myVideo, frame);
    end
    
    % update external forces vectors
    fu(BC.fixed_u) = Solution.fuE;
    fp(BC.fixed_p) = Solution.fpE;
	ff(BC.fixed_uf) = Solution.ffE;
    
    % post processing: compute stress/flux, export VTK file
    if plot2vtk
        PostProcessing(Solution, Material, MeshU, MeshP, MeshN, Control, BC, config_name, vtk_dir);
    end

    % store pressure over time
    Plot.p_time(Control.step,:) = Solution.p(Control.plotp, 1);
    % store solid displacement over time
    Plot.u_time(Control.step,:) = Solution.u(Control.plotu, 1);
    % store solid velocity over time
    Plot.udot_time(Control.step,:) = Solution.udot(Control.plotu, 1);
    % store fluid displacement over time
    Plot.uf_time(Control.step,:) = Solution.uf(Control.plotu, 1);
    % store fluid velocity over time
    Plot.ufdot_time(Control.step,:) = Solution.ufdot(Control.plotu, 1);
    % store solid acceleration over time
    Plot.u2dot_time(Control.step,:) = Solution.u2dot(Control.plotu, 1);
    % store fluid acceleration over time
    Plot.uf2dot_time(Control.step,:) = Solution.uf2dot(Control.plotu, 1);

    % store variables over time at fixed row
    Plot.u_synthetic(Control.step,:) = Solution.u(Control.ploturow);
    Plot.udot_synthetic(Control.step,:) = Solution.udot(Control.ploturow);
    Plot.p_synthetic(Control.step,:) = Solution.p(Control.plotprow);
    Plot.uf_synthetic(Control.step,:) = Solution.uf(Control.ploturow);
    Plot.ufdot_synthetic(Control.step,:) = Solution.ufdot(Control.ploturow);

    % store variables over space
    Plot.p_space = Solution.p(Control.plotp);
    Plot.u_space = Solution.u(Control.plotu);
    Plot.udot_space = Solution.udot(Control.plotu);
    Plot.uf_space = Solution.uf(Control.plotu);
    Plot.ufdot_space = Solution.ufdot(Control.plotu);
    Plot.u2dot_space = Solution.u2dot(Control.plotu);
    Plot.uf2dot_space = Solution.uf2dot(Control.plotu);
    
    % update variables
    Iteration.u_old = Solution.u; % solid displacement
    Iteration.udot_old = Solution.udot; % solid velocity
    Iteration.u2dot_old = Solution.u2dot; % solid acceleration
    Iteration.p_old = Solution.p; % fluid pressure
    Iteration.pdot_old = Solution.pdot; % fluid pressure gradient
    Iteration.uf_old = Solution.uf; % fluid displacement
    Iteration.ufdot_old = Solution.ufdot; % fluid velocity
    Iteration.uf2dot_old = Solution.uf2dot; % fluid acceleration
    Iteration.fu_old = fu; % load vector
    Iteration.ff_old = ff; % load vector
end
 
% store variables at fixed row
Plot.urow = Solution.u(Control.ploturow);
Plot.udotrow = Solution.udot(Control.ploturow);
Plot.prow = Solution.p(Control.plotprow);
Plot.ufrow = Solution.uf(Control.ploturow);
Plot.ufdotrow = Solution.ufdot(Control.ploturow);

% close video file
if saveVideo_on
    close(myVideo)
end
