% Biot poroelasticity model
% Dynamic case
% December 2022
% ------------------------------------------------------------------------

%% Model name and type
disp([num2str(toc),': Model: Biot dynamic case']);

%% Assemble system matrices
disp([num2str(toc),': Assembling System Matrices...']);

[Kuu, Kup, Kpp, M, Mhat, S] = ComputeMatricesDynamic_Biot(Material, MeshU, MeshP, QuadU, QuadP);

%% Assemble system load vectors
[fu,fp,fn] = ComputeLoads(BC, MeshU, MeshP, MeshN, Control, QuadU, QuadP);

%% Solve eigenproblem
if Control.freqDomain
    disp([num2str(toc),': Solving Uncoupled Eigenproblems...']);
    [phi_u, omega2_u, phi_p, omega2_p] = SolveEig_Dynamic(Kuu, Kup, Kpp, M, Mhat, MeshU, MeshP, BC, Control);
else
    phi_u = [];
    phi_p = [];
end

%% Initialize iteration variables
[Iteration, Plot] = initVariables(phi_u, phi_p, MeshU, MeshP, MeshN, Material, Control, BC);

%% Initial condition file
if plot2vtk
    % time domain solution
    Solution.u = Iteration.u_old;
    Solution.udot = Iteration.udot_old;
    Solution.u2dot = Iteration.u2dot_old;
    Solution.p = Iteration.p_old;

    % save vtk file
    PostProcessing(Solution, Material, MeshU, MeshP, MeshN, Control, BC, config_name, vtk_dir);
    % update time step
    Control.step = 1;
end

%% Solve system
for t = 0:Control.dt:Control.tend
    fprintf('\n Step %d, t = %d \n', Control.step, t);
    Control.t = t;

    % linear solver
    [Solution] = SolverDynamic_Biot(Kuu, Kup, Kpp, M, Mhat, S, fu, fp, BC, Control, Iteration);

%     figure(1);
%     subplot(1,2,1);
%     plot(MeshU.coords, Solution.u, 'm', 'LineWidth', 1.5);
%     title('Displacement');
%     xlim([0 max(MeshU.coords)]);
%     subplot(1,2,2);
%     plot(MeshU.coords, Solution.udot, 'b', 'LineWidth', 1.5);
%     title('Velocity');
%     xlim([0 max(MeshU.coords)]);
%     pause(0.001);

    % solution in the frequency domain
    if Control.freqDomain
        [SolutionFreq] = SolverDynamicFreq(phi_u, omega2_u, phi_p, omega2_p, Kuu, Kup, Kpp, M, Mhat, S, fu, fp, BC, Control, Iteration);
    end

    % update external forces vectors
    fu(BC.fixed_u) = Solution.fE;
    fp(BC.fixed_p) = Solution.qE;

    % post processing: compute stress/flux, export VTK file
    if plot2vtk
        PostProcessing(Solution, Material, MeshU, MeshP, MeshN, Control, BC, config_name, vtk_dir);
    end

    % store variables over time
    if Control.step < length(Plot.time)
        % time domain
        % plot pressure vs time
        Plot.p_time(Control.step+1,:) = Solution.p(Control.plotp, 1);
        % plot displacement vs time
        Plot.u_time(Control.step+1,:) = Solution.u(Control.plotu, 1);
        % plot velocity vs time
        Plot.udot_time(Control.step+1,:) = Solution.udot(Control.plotu, 1);

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

    % update variables - time domain
    Iteration.u_old = Solution.u; % solid displacement
    Iteration.udot_old = Solution.udot; % solid velocity
    Iteration.p_old = Solution.p; % fluid pressure
    Iteration.u2dot_old = Solution.u2dot; % solid acceleration
    Iteration.pdot_old = Solution.pdot; % fluid pressure gradient
    Iteration.fu_old = fu; % load vector
    Iteration.fp_old = fp; % flux vector

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
end
