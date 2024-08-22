% de la Cruz and Spanos poroelasticity model
% Dynamic case - Mode superposition
% December 2022
% ------------------------------------------------------------------------

%% Model name and type
disp([num2str(toc),': Model: Spanos u-p-n dynamic case, mode superposition']);

%% Assemble system matrices
disp([num2str(toc),': Assembling System Matrices...']);

[Muu, Mpu, Mnu, Kuu, Kup, Kpp, Kpu, S, Kpn, Knn, Knu, Knp, Kun] = ComputeMatricesDyn2_Spanos_UPN(Material, MeshU, MeshP, MeshN, QuadU);

%% Solve eigenproblem
disp([num2str(toc),': Solving Uncoupled Eigenproblems...']);
[phi_u, omega2_u, phi_p, omega2_p, phi_n, omega2_n] = EigenDyn_UPN(Muu, Mpu, Mnu, Kuu, Kup, Kpp, Knp, MeshU, MeshP, MeshN, BC, Control);

%% Initialize iteration variables
[Iteration, Plot] = initVariables_Freq(phi_u, phi_p, phi_n, MeshU, MeshP, MeshN, Material, Control, BC);

%% Initial condition file
if plot2vtk
    % time domain solution
    Solution.u = Iteration.u_old;
    Solution.udot = Iteration.udot_old;
    Solution.u2dot = Iteration.u2dot_old;
    Solution.p = Iteration.p_old;
    Solution.n = Iteration.n_old;

    % save vtk file
    PostProcessing(Solution, Material, MeshU, MeshP, MeshN, Control, BC, config_name, vtk_dir);

    % update time step
    Control.step = 1;
end

%% Solve system
for t = 1:length(Plot.time)
    % current time
    Control.t = Plot.time(t);
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
    [fu,fp,fn] = ComputeLoads(BC, MeshU, MeshP, MeshN, Control, QuadU, QuadP, Material);

    % linear solver
    [Solution] = SolverDyn_UPN(Muu, Mpu, Mnu, Kuu, Kup, Kpp, Kpu, S, Kpn, Knn, Knu, Knp, Kun, fu, fp, fn, BC, Control, Iteration);

    % solution using mode superposition
    [SolutionFreq] = SolverFreqDyn_UPN(phi_u, omega2_u, phi_p, omega2_p, phi_n, omega2_n, Kuu, Kup, Kpu, Kpp, S, Kpn, Knu, Knp, Knn, Muu, Mpu, Mnu, fu, fp, BC, Control, Iteration);

    % update external forces vectors
    fu(BC.fixed_u) = Solution.fE;
    fp(BC.fixed_p) = Solution.qE;
    fn(BC.fixed_n) = Solution.fnE;

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
        % plot porosity vs time
        Plot.n_time(Control.step+1,:) = Solution.n(Control.plotp, 1);
        % plot solid acceleration vs time
        Plot.u2dot_time(Control.step+1,:) = Solution.u2dot(Control.plotu, 1);
        
        % mode superposition
        % plot pressure vs time
        Plot.pF(Control.step+1,:) = SolutionFreq.pF(Control.plotp, 1);
        % plot displacement vs time
        Plot.uF(Control.step+1,:) = SolutionFreq.uF(Control.plotu, 1);
        % plot velocity vs time
        Plot.uFdot(Control.step+1,:) = SolutionFreq.uFdot(Control.plotu, 1);
        % plot porosity vs time
        Plot.nF(Control.step+1,:) = SolutionFreq.nF(Control.plotp, 1);

        % synthetics
        Plot.u_synthetic(Control.step+1,:) = Solution.u(Control.ploturow);
        Plot.udot_synthetic(Control.step+1,:) = Solution.udot(Control.ploturow);
        Plot.p_synthetic(Control.step+1,:) = Solution.p(Control.plotprow);
        Plot.n_synthetic(Control.step+1,:) = Solution.n(Control.plotprow);
    end

    % store variables over space
    Plot.p_space = Solution.p(Control.plotp);
    Plot.u_space = Solution.u(Control.plotu);
    Plot.udot_space = Solution.udot(Control.plotu);
    Plot.u2dot_space = Solution.u2dot(Control.plotu);

    % update variables - time domain
    Iteration.u_old = Solution.u; % solid displacement
    Iteration.udot_old = Solution.udot; % solid velocity
    Iteration.p_old = Solution.p; % fluid pressure
    Iteration.n_old = Solution.n; % medium porosity
    Iteration.u2dot_old = Solution.u2dot; % solid acceleration
    Iteration.pdot_old = Solution.pdot; % fluid pressure gradient
    Iteration.ndot_old = Solution.ndot; % porosity variation
    Iteration.fu_old = fu; % load vector
    Iteration.fp_old = fp; % flux vector
    Iteration.fn_old = fn; % flux vector

    % update variables - mode superposition
    Iteration.xuF_old = SolutionFreq.xuF;
    Iteration.xuFdot_old = SolutionFreq.xuFdot;
    Iteration.xpF_old = SolutionFreq.xpF;
    Iteration.xuF2dot_old = SolutionFreq.xuF2dot;
    Iteration.xpFdot_old = SolutionFreq.xpFdot;
    Iteration.xnF_old = SolutionFreq.xnF;
    Iteration.xnFdot_old = SolutionFreq.xnFdot;

    Iteration.uF_old = SolutionFreq.uF; % solid displacement
    Iteration.uFdot_old = SolutionFreq.uFdot; % solid velocity
    Iteration.pF_old = SolutionFreq.pF; % fluid pressure
    Iteration.uF2dot_old = SolutionFreq.uF2dot; % solid acceleration
    Iteration.pFdot_old = SolutionFreq.pFdot; % fluid pressure gradient
    Iteration.nF_old = SolutionFreq.nF; % porosity
    Iteration.nFdot_old = SolutionFreq.nFdot; % porosity gradient

    % update time step
    Control.step = Control.step + 1;
end

% plot variables in length for fixed coordinate
Plot.urow = Solution.u(Control.ploturow);
Plot.udotrow = Solution.udot(Control.ploturow);
Plot.prow = Solution.p(Control.plotprow);
