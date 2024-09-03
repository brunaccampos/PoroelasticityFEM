% Biot poroelasticity model
% Dynamic case
% December 2022
% ------------------------------------------------------------------------

%% Model name and type
disp([num2str(toc),': Model: Biot u-p-w dynamic case']);

%% Assemble system matrices
disp([num2str(toc),': Assembling System Matrices...']);
[Mss, Msf, Kss, Ksp, Mfs, Mff, Kfp, Kff, Cps, Kpf, Cpp, Cfs] = ComputeMatricesDyn7_Biot_UPW(Material, MeshU, MeshP, QuadU, QuadP);

%% Initial condition file
if plot2vtk
    % time domain solution
    Solution.u = Iteration.u_old;
    Solution.udot = Iteration.udot_old;
    Solution.u2dot = Iteration.u2dot_old;
    Solution.p = Iteration.p_old;
    Solution.w = Iteration.w_old;
    Solution.wdot = Iteration.wdot_old;
    
    % save vtk file
    PostProcessing(Solution, Material, MeshU, MeshP, MeshN, Control, BC, config_name, vtk_dir);
end
    % update time step
    Control.step = 1;

% initialize video file
if saveVideo_on
    myVideo = VideoWriter('myVideoFile'); % open video file
    myVideo.FrameRate = 20;
    open(myVideo)
end

% print time and step
fprintf('Step 0, t = 0');
msg_len = numel('Step 0, t = 0');

%% Solve system
for t = 1:length(Plot.time)
    % current time
    Control.t = Plot.time(t);
    % print current time and step
    fprintf(repmat('\b',1,msg_len));
    msg_len = fprintf('Step %d, t = %.4g',Control.step', Control.t);
    
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
    [fu, fp, ff] = ComputeLoads_UPW(BC, MeshU, MeshP, Control, Material, QuadU, QuadP);
    
    % analytical solution for 1D case
    if Control.plotansol
        [p_an, u_an, w_an] = feval(Control.ansol_type, Control, Material, MeshU, MeshP, BC);
       
        % store variables over length
        Plot.pan_space = p_an;
        Plot.uan_space = u_an;
        Plot.wan_space = w_an;
        
        % store variables over time
        if Control.step < length(Plot.time)
            % store variables over time
            Plot.uan_time(Control.step+1,:) = u_an(Control.plotu, 1);
            Plot.pan_time(Control.step+1,:) = p_an(Control.plotp, 1);
            Plot.wan_time(Control.step+1,:) = w_an(Control.plotu, 1);
        end
    end
    
    % linear solver
    [Solution] = SolverDyn_UPW(Mss, Msf, Kss, Ksp, Mfs, Mff, Kfp, Kff, Cps, Kpf, Cpp, Cfs, fu, fp, ff, BC, Control, Iteration);

    % update external forces vectors
    fu(BC.fixed_u) = Solution.fuE;
	ff(BC.fixed_w) = Solution.ffE;
    fp(BC.fixed_p) = Solution.fpE;
    
    % post processing: compute stress/flux, export VTK file
    if plot2vtk
        PostProcessing(Solution, Material, MeshU, MeshP, MeshN, Control, BC, config_name, vtk_dir);
    end

    % store variables over time
    if Control.step < length(Plot.time)
        % time domain
        % plot pressure vs time
        Plot.p_time(Control.step+1,:) = Solution.p(Control.plotp, 1);
        % plot solid displacement vs time
        Plot.u_time(Control.step+1,:) = Solution.u(Control.plotu, 1);
        % plot solid velocity vs time
        Plot.udot_time(Control.step+1,:) = Solution.udot(Control.plotu, 1);
        % plot relative fluid velocity vs time
        Plot.w_time(Control.step+1,:) = Solution.w(Control.plotu, 1);
        % plot solid acceleration vs time
        Plot.u2dot_time(Control.step+1,:) = Solution.u2dot(Control.plotu, 1);
        % plot relative fluid acceleration vs time
        Plot.wdot_time(Control.step+1,:) = Solution.wdot(Control.plotu, 1);
        
        % synthetics
        Plot.u_synthetic(Control.step+1,:) = Solution.u(Control.ploturow);
        Plot.udot_synthetic(Control.step+1,:) = Solution.udot(Control.ploturow);
        Plot.p_synthetic(Control.step+1,:) = Solution.p(Control.plotprow);
        Plot.w_synthetic(Control.step+1,:) = Solution.w(Control.ploturow);
        Plot.wdot_synthetic(Control.step+1,:) = Solution.wdot(Control.ploturow);
    end

    % store variables over space
    Plot.p_space = Solution.p(Control.plotp);
    Plot.u_space = Solution.u(Control.plotu);
    Plot.udot_space = Solution.udot(Control.plotu);
    Plot.u2dot_space = Solution.u2dot(Control.plotu);
    Plot.w_space = Solution.w(Control.plotu);
    Plot.wdot_space = Solution.wdot(Control.plotu);
    
    % update variables - time domain
    Iteration.u_old = Solution.u; % solid displacement
    Iteration.udot_old = Solution.udot; % solid velocity
    Iteration.u2dot_old = Solution.u2dot; % solid acceleration
   
    Iteration.p_old = Solution.p; % fluid pressure
    Iteration.pdot_old = Solution.pdot; % fluid pressure gradient
    
    Iteration.w_old = Solution.w; % relative fluid velocity
    Iteration.wdot_old = Solution.wdot; % relative fluid acceleration
    
    Iteration.fu_old = fu; % load vector
    Iteration.ff_old = ff; % load vector
    Iteration.fp_old = fp; % load vector
    
    % update time step
    Control.step = Control.step + 1;
end
 
% plot variables in length for fixed coordinate
Plot.urow = Solution.u(Control.ploturow);
Plot.udotrow = Solution.udot(Control.ploturow);
Plot.prow = Solution.p(Control.plotprow);
Plot.wrow = Solution.w(Control.ploturow);

% close video file
if saveVideo_on
    close(myVideo)
end