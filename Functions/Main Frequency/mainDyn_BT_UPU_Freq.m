% Biot poroelasticity model
% Dynamic case - Mode superposition
% December 2022
% ------------------------------------------------------------------------

%% Model name and type
disp([num2str(toc),': Model: BT u-p-U dynamic, mode superposition']);

%% Assemble system matrices
disp([num2str(toc),': Assembling System Matrices...']);
[Kss, Ksp, Mss, Csf, Css, Kpf, Kps, Kpp, Kfp, Mff, Cff, Cfs, Msf, Mfs] = ComputeMatricesDyn_BT_UPU(Material, MeshU, MeshP, QuadU, QuadP);

%% Solve eigenproblem
disp([num2str(toc),': Solving Uncoupled Eigenproblems...']);
[phi_u, omega2_u, phi_p, omega2_p, phi_uf, omega2_uf] = EigenDyn_UPU(Kss, Ksp, Mss, Kpf, Kps, Kpp, Kfp, Mff, Msf, Mfs, MeshU, MeshP, BC, Control);

%% Initialize iteration variables
[Iteration, Plot] = initVariables_Freq(phi_u, phi_p, [], MeshU, MeshP, MeshN, Material, Control, BC);

%% Initial condition file
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
    % update time step
    Control.step = 1;
end

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
    [fu, fp, ff] = ComputeLoads_UPU(BC, MeshU, MeshP, Control, Material, QuadU, QuadP);
    
    % analytical solution for 1D case
    if Control.plotansol
        [p_an, u_an, uf_an] = feval(Control.ansol_type, Control, Material, MeshU, MeshP, BC);
       
        % store variables over length
        Plot.pan_space = p_an;
        Plot.uan_space = u_an;
        Plot.ufan_space = uf_an;
        
        % store variables over time
        if Control.step < length(Plot.time)
            % store variables over time
            Plot.uan_time(Control.step+1,:) = u_an(Control.plotu, 1);
            Plot.pan_time(Control.step+1,:) = p_an(Control.plotp, 1);
            Plot.ufan_time(Control.step+1,:) = uf_an(Control.plotu, 1);
        end
    end
    
    % linear solver
    [Solution] = SolverDyn_UPU(Kss, Ksp, Mss, Csf, Css, Kpf, Kps, Kpp, Kfp, Mff, Cff, Cfs, Msf, Mfs, fu, fp, ff, BC, Control, Iteration);
    
    % solution using mode superposition
    [SolutionFreq] = SolverFreqDyn_UPU(phi_u, omega2_u, phi_p, omega2_p, phi_uf, omega_uf, Kss, Ksp, Mss, Csf, Css, Kpf, Kps, Kpp, Kfp, Mff, Cff, Cfs, Msf, Mfs, fu, ff, BC, Control, Iteration);

    % update external forces vectors
    fu(BC.fixed_u) = Solution.fuE;
    fp(BC.fixed_p) = Solution.fpE;
	ff(BC.fixed_uf) = Solution.ffE;
    
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
        % plot fluid displacement vs time
        Plot.uf_time(Control.step+1,:) = Solution.uf(Control.plotu, 1);
        % plot fluid velocity vs time
        Plot.ufdot_time(Control.step+1,:) = Solution.ufdot(Control.plotu, 1);
        % plot solid acceleration vs time
        Plot.u2dot_time(Control.step+1,:) = Solution.u2dot(Control.plotu, 1);
        % plot fluid acceleration vs time
        Plot.uf2dot_time(Control.step+1,:) = Solution.uf2dot(Control.plotu, 1);
        
        % mode superposition
        % plot pressure vs time
        Plot.pF(Control.step+1,:) = SolutionFreq.pF(Control.plotp, 1);
        % plot displacement vs time
        Plot.uF(Control.step+1,:) = SolutionFreq.uF(Control.plotu, 1);
        % plot velocity vs time
        Plot.uFdot(Control.step+1,:) = SolutionFreq.uFdot(Control.plotu, 1);

        % synthetics
        Plot.u_synthetic(Control.step+1,:) = Solution.u(Control.ploturow);
        Plot.udot_synthetic(Control.step+1,:) = Solution.udot(Control.ploturow);
        Plot.p_synthetic(Control.step+1,:) = Solution.p(Control.plotprow);
        Plot.uf_synthetic(Control.step+1,:) = Solution.uf(Control.ploturow);
        Plot.ufdot_synthetic(Control.step+1,:) = Solution.ufdot(Control.ploturow);
    end

    % store variables over space
    Plot.p_space = Solution.p(Control.plotp);
    Plot.u_space = Solution.u(Control.plotu);
    Plot.udot_space = Solution.udot(Control.plotu);
    Plot.uf_space = Solution.uf(Control.plotu);
    Plot.ufdot_space = Solution.ufdot(Control.plotu);
    Plot.u2dot_space = Solution.u2dot(Control.plotu);
    Plot.uf2dot_space = Solution.uf2dot(Control.plotu);
    
    % update variables - time domain
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

    % update variables - mode superposition
    Iteration.xuF_old = SolutionFreq.xuF;
    Iteration.xuFdot_old = SolutionFreq.xuFdot;
    Iteration.xpF_old = SolutionFreq.xpF;
    Iteration.xuF2dot_old = SolutionFreq.xuF2dot;
    Iteration.xpFdot_old = SolutionFreq.xpFdot;

    Iteration.uF_old = SolutionFreq.uF; % solid displacement
    Iteration.uFdot_old = SolutionFreq.uFdot; % solid velocity
    Iteration.pF_old = SolutionFreq.pF; % fluid pressure
    Iteration.uF2dot_old = SolutionFreq.uF2dot; % solid acceleration
    Iteration.pFdot_old = SolutionFreq.pFdot; % fluid pressure gradient

    % update time step
    Control.step = Control.step + 1;
end
 
% plot variables in length for fixed coordinate
Plot.urow = Solution.u(Control.ploturow);
Plot.udotrow = Solution.udot(Control.ploturow);
Plot.prow = Solution.p(Control.plotprow);
Plot.ufrow = Solution.uf(Control.ploturow);
Plot.ufdotrow = Solution.ufdot(Control.ploturow);

% close video file
if saveVideo_on
    close(myVideo)
end