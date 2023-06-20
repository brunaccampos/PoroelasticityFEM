% de la Cruz and Spanos poroelasticity model
% Transient case
% December 2022
% ------------------------------------------------------------------------

%% Model name and type
disp([num2str(toc),': Model: Spanos u-p-n transient case']);

%% Assemble system matrices
disp([num2str(toc),': Assembling System Matrices...']);

[Kuu, Kup, Kpp, Kpu, S, Kpn, Knn, Knu, Knp, Kun] = ComputeMatricesTr2_Spanos_UPN(Material, MeshU, MeshP, MeshN, QuadU, QuadP);

%% Solve eigenproblem
if Control.freqDomain
    disp([num2str(toc),': Solving Uncoupled Eigenproblems...']);
    [phi_u, omega2_u, phi_p, omega2_p, phi_n, omega2_n] = SolveEigTransient_Spanos(Kuu, Kup, Kpp, Knp, Kpu, S, Kpn, Knu, Knn, MeshU, MeshP, MeshN, BC, Control);
else
    phi_u = [];
    phi_p = [];
    phi_n = [];
end

%% Initialize iteration variables
[Iteration, Plot] = initVariables(phi_u, phi_p, phi_n, MeshU, MeshP, MeshN, Material, Control, BC);

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

% initialize time variable
Control.t = 0;

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
    [fu,fp,fn] = ComputeLoads(BC, MeshU, MeshP, MeshN, Control, QuadU, QuadP, Material);

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
        if Control.step < length(Plot.time)
            Plot.uan_time(Control.step+1,:) = u_an(Control.plotu, 1);
            Plot.pan_time(Control.step+1,:) = p_an(Control.plotp, 1);
        end
    end

    % linear solver
    [Solution] = SolverTr_UPN(Kuu, Kup, Kpp, Kpu, S, Kpn, Knn, Knu, Knp, Kun, fu, fp, fn, BC, Control, Iteration);
    
    % solution in the frequency domain
    if Control.freqDomain
        [SolutionFreq] = SolverTransientFreq_Spanos(phi_u, omega2_u, phi_p, omega2_p, phi_n, omega2_n, Kuu, Kup, Kpu, Kpp, S, Kpn, Knu, Knp, Knn, fu, fp, BC, Control, Iteration);
    end
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
        
        % frequency domain
        if Control.freqDomain
            % plot pressure vs time
            Plot.pF(Control.step+1,:) = SolutionFreq.pF(Control.plotp, 1);
            % plot displacement vs time
            Plot.uF(Control.step+1,:) = SolutionFreq.uF(Control.plotu, 1);
            % plot velocity vs time
            Plot.uFdot(Control.step+1,:) = SolutionFreq.uFdot(Control.plotu, 1);
            % plot porosity vs time
            Plot.nF(Control.step+1,:) = SolutionFreq.nFdot(Control.plotp, 1);
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
    Iteration.n_old = Solution.n; % medium porosity
    Iteration.fu_old = fu; % load vector
    Iteration.fp_old = fp; % flux vector  
    Iteration.fn_old = fn; % flux vector  

    % update variables - frequency domain
    if Control.freqDomain
        Iteration.xuF_old = SolutionFreq.xuF;
        Iteration.xuFdot_old = SolutionFreq.xuFdot;
        Iteration.xpF_old = SolutionFreq.xpF;
        Iteration.xnF_old = SolutionFreq.xnF;
        Iteration.xnFdot_old = SolutionFreq.xnFdot;

        Iteration.uF_old = SolutionFreq.uF; % solid displacement
        Iteration.uFdot_old = SolutionFreq.uFdot; % solid velocity
        Iteration.pF_old = SolutionFreq.pF; % fluid pressure
        Iteration.nF_old = SolutionFreq.nF; % porosity
        Iteration.nFdot_old = SolutionFreq.nFdot; % porosity change
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
    Plot.nrow = Solution.n(Control.plotprow);
end
