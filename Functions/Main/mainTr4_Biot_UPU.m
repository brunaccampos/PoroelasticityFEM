% Biot poroelasticity model
% Transient case
% September 2023
% ------------------------------------------------------------------------

%% Model name and type
disp([num2str(toc),': Model: Biot u-p-U transient case']);

%% Assemble system matrices
disp([num2str(toc),': Assembling System Matrices...']);

[Kss, Ksp, Csf, Css, Kpf, Kps, Kpp, Kfp, Cff, Cfs] = ComputeMatricesTr4_Biot_UPU(Material, MeshU, MeshP, QuadU, QuadP);

%% Solve eigenproblem
if Control.freqDomain
    disp([num2str(toc),': Solving Uncoupled Eigenproblems...']);
%     [phi_u, omega2_u, phi_p, omega2_p, phi_uf, omega2_uf] = EigenTr_UPU(Kss, Ksp, Kpf, Kps, Kpp, Kfp, MeshU, MeshP, BC, Control);
else
    phi_u = [];
    phi_p = [];
end

%% Initialize iteration variables
[Iteration, Plot] = initVariables(phi_u, phi_p, [], MeshU, MeshP, MeshN, Material, Control, BC);

%% Initial condition file
if plot2vtk
    % time domain solution
    Solution.u = Iteration.u_old;
    Solution.udot = Iteration.udot_old;
    Solution.p = Iteration.p_old;
    Solution.uf = Iteration.uf_old;
    Solution.ufdot = Iteration.ufdot_old;
    
    % save vtk file
    PostProcessing(Solution, Material, MeshU, MeshP, MeshN, Control, BC, config_name, vtk_dir);
   
    % update time step
    Control.step = 1;
end

% initialize time variable
Control.t = 0;

% initialize video file
if saveVideo_on
    myVideo = VideoWriter('myVideoFile'); %open video file
    myVideo.FrameRate = 20;
    open(myVideo)
end

%% Solve system
while Control.t < Control.tend
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
        if Control.uncoupled
            p_an = Control.p_an(Control.t);
            u_an = Control.u_an(Control.t);
        else
            if any(Material.Minv)
                [p_an, u_an] = getAnalyticResult_Comp(Material, MeshU, MeshP, BC, Control);
            else % 1/M = 0
                [~, p_an, u_an] = getAnalyticResult_Incomp(Material, MeshU, MeshP, BC, Control);
            end
        end
       
        % store variables over length
        Plot.pan_space = p_an;
        Plot.uan_space = u_an;

        % store variables over time
        Plot.uan_time(Control.step,:) = u_an(Control.plotu, 1);
        Plot.pan_time(Control.step,:) = p_an(Control.plotp, 1);
    end
    
    % linear solver
    [Solution] = SolverTr_UPU(Kss, Ksp, Csf, Css, Kpf, Kps, Kpp, Kfp, Cff, Cfs, fu, fp, ff, BC, Control, Iteration);

    % plot solution over time
    figure(1);
    % solid displacement
    subplot(2,3,1);
    grid on
    plot(MeshU.coords, Solution.u, 'm', 'LineWidth', 1.5);
    title('Solid displacement');
    hold on
    plot(MeshU.coords, Plot.uan_space, 'k--', 'LineWidth', 1.5);
    hold off
    % fluid pressure
    subplot(2,3,2);
    grid on
    plot(MeshP.coords, Solution.p, 'g', 'LineWidth', 1.5);
    title('Pressure');
    hold on
    plot(MeshP.coords, Plot.pan_space, 'k--', 'LineWidth', 1.5);
    hold off
    % fluid displacement
    subplot(2,3,3);
    grid on
    plot(MeshU.coords, Solution.uf, 'm', 'LineWidth', 1.5);
    title('Fluid displacement');
    hold on
    plot(MeshU.coords, Plot.uan_space, 'k--', 'LineWidth', 1.5);
    hold off
    pause(0.001);
    if saveVideo_on
        frame = getframe(gcf); %get frame
        writeVideo(myVideo, frame);
    end
    
    % solution in the frequency domain
    if Control.freqDomain
%         [SolutionFreq] = SolverFreqTr_UPU(phi_u, omega2_u, phi_p, omega2_p, phi_uf, omega_uf, Kss, Ksp, Mss, Csf, Css, Kpf, Kps, Kpp, Kfp, Mff, Cff, Cfs, Msf, Mfs, fu, ff, BC, Control, Iteration);
    end

    % update external forces vectors
    fu(BC.fixed_u) = Solution.fuE;
    fp(BC.fixed_p) = Solution.fpE;
	ff(BC.fixed_u) = Solution.ffE;
    
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
        
        % frequency domain
        if Control.freqDomain
            % plot pressure vs time
            Plot.pF(Control.step+1,:) = SolutionFreq.pF(Control.plotp, 1);
            % plot displacement vs time
            Plot.uF(Control.step+1,:) = SolutionFreq.uF(Control.plotu, 1);
            % plot velocity vs time
            Plot.uFdot(Control.step+1,:) = SolutionFreq.uFdot(Control.plotu, 1);
        end

        % synthetics
        if Control.fixedDepthPlotON
            Plot.u_synthetic(Control.step+1,:) = Solution.u(Control.ploturow);
            Plot.udot_synthetic(Control.step+1,:) = Solution.udot(Control.ploturow);
            Plot.p_synthetic(Control.step+1,:) = Solution.p(Control.plotprow);
            Plot.uf_synthetic(Control.step+1,:) = Solution.uf(Control.ploturow);
            Plot.ufdot_synthetic(Control.step+1,:) = Solution.ufdot(Control.ploturow);
        end
    end

    % store variables over space
    Plot.p_space = Solution.p(Control.plotp);
    Plot.u_space = Solution.u(Control.plotu);
    Plot.udot_space = Solution.udot(Control.plotu);
    Plot.uf_space = Solution.uf(Control.plotu);
    Plot.ufdot_space = Solution.ufdot(Control.plotu);
    
    % update variables - time domain
    Iteration.u_old = Solution.u; % solid displacement
    Iteration.udot_old = Solution.udot; % solid velocity
   
    Iteration.p_old = Solution.p; % fluid pressure
    Iteration.pdot_old = Solution.pdot; % fluid pressure gradient
    
    Iteration.uf_old = Solution.uf; % fluid displacement
    Iteration.ufdot_old = Solution.ufdot; % fluid velocity
    
    Iteration.fu_old = fu; % load vector
    Iteration.fp_old = fp; % load vector
    Iteration.ff_old = ff; % load vector

    % update variables - frequency domain
    if Control.freqDomain
        Iteration.xuF_old = SolutionFreq.xuF;
        Iteration.xuFdot_old = SolutionFreq.xuFdot;
        Iteration.xpF_old = SolutionFreq.xpF;
        Iteration.xpFdot_old = SolutionFreq.xpFdot;

        Iteration.uF_old = SolutionFreq.uF; % solid displacement
        Iteration.uFdot_old = SolutionFreq.uFdot; % solid velocity
        Iteration.pF_old = SolutionFreq.pF; % fluid pressure
        Iteration.pFdot_old = SolutionFreq.pFdot; % fluid pressure gradient
    end

    % update time and step
    Control.step = Control.step + 1;
    Control.t = Control.t + Control.dtc;
end
 
% plot variables in x for fixed y (2D case)
if Control.fixedDepthPlotON
    Plot.urow = Solution.u(Control.ploturow);
    Plot.udotrow = Solution.udot(Control.ploturow);
    Plot.prow = Solution.p(Control.plotprow);
    Plot.ufrow = Solution.uf(Control.ploturow);
    Plot.ufdotrow = Solution.ufdot(Control.ploturow);
end

% close video file
if saveVideo_on
    close(myVideo)
end
