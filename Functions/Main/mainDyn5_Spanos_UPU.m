% Biot poroelasticity model
% Dynamic case
% December 2022
% ------------------------------------------------------------------------

%% Model name and type
disp([num2str(toc),': Model: Spanos u-p-U dynamic case']);

%% Assemble system matrices
disp([num2str(toc),': Assembling System Matrices...']);

[Kss, Ksp, Mss, Csf, Css, Kpf, Kps, Kpp, Kfp, Mff, Cff, Cfs] = ComputeMatricesDyn5_Spanos_UPU(Material, MeshU, MeshP, QuadU, QuadP);

%% Solve eigenproblem
if Control.freqDomain
    disp([num2str(toc),': Solving Uncoupled Eigenproblems...']);
%     [phi_u, omega2_u, phi_p, omega2_p] = SolveEigDynamic_Biot(Kuu, Kup, Kpp, Muu, Mpu, MeshU, MeshP, BC, Control);
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

% initialize time variable
Control.t = 0;

%% Solve system
while Control.t < Control.tend
    fprintf('\n Step %d, t = %d \n', Control.step, Control.t);

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
    [fu, ff] = ComputeLoads_UPU(BC, MeshU, MeshP, Control, Material, QuadU, QuadP);
    
    % linear solver
    [Solution] = SolverDyn_UPU(Kss, Ksp, Mss, Csf, Css, Kpf, Kps, Kpp, Kfp, Mff, Cff, Cfs, fu, ff, BC, Control, Iteration);

    % solution in the frequency domain
    if Control.freqDomain
%         [SolutionFreq] = SolverDynamicFreq_Biot(phi_u, omega2_u, phi_p, omega2_p, Kuu, Kup, Kpp, Muu, Mpu, S, fu, fp, BC, Control, Iteration);
    end

    % update external forces vectors
    fu(BC.fixed_u) = Solution.fuE;
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
        % plot solid acceleration vs time
        Plot.u2dot_time(Control.step+1,:) = Solution.u2dot(Control.plotu, 1);
        % plot fluid acceleration vs time
        Plot.uf2dot_time(Control.step+1,:) = Solution.uf2dot(Control.plotu, 1);
        
        % frequency domain
        if Control.freqDomain
            % plot pressure vs time
            Plot.pF(Control.step+1,:) = SolutionFreq.pF(Control.plotp, 1);
            % plot displacement vs time
            Plot.uF(Control.step+1,:) = SolutionFreq.uF(Control.plotu, 1);
            % plot velocity vs time
            Plot.uFdot(Control.step+1,:) = SolutionFreq.uFdot(Control.plotu, 1);
        end
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

    % update variables - frequency domain
    if Control.freqDomain
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
    end

    % update time and step
    Control.step = Control.step + 1;
    Control.t = Control.t + Control.dtc;
end
 
% plot variables in x for fixed y (2D case)
if isfield(Control, 'depthplot')
    Plot.urow = Solution.u(Control.ploturow);
    Plot.udotrow = Solution.udot(Control.ploturow);
    Plot.prow = Solution.p(Control.plotprow);
    Plot.ufrow = Solution.uf(Control.ploturow);
    Plot.ufdotrow = Solution.ufdot(Control.ploturow);
end
