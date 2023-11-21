% Biot poroelasticity model
% Transient case
% December 2022
% ------------------------------------------------------------------------

%% Model name and type
disp([num2str(toc),': Model: Biot u-p transient case']);

%% Assemble system matrices
disp([num2str(toc),': Assembling System Matrices...']);

[Kuu, Kup, Kpp, Kpu, S] = ComputeMatricesTr1_Biot_UP(Material, MeshU, MeshP, QuadU, QuadP);

%% Solve eigenproblem
if Control.freqDomain
    disp([num2str(toc),': Solving Uncoupled Eigenproblems...']);
    [phi_u, omega2_u, phi_p, omega2_p] = EigenTr_UP(Kuu, Kup, Kpp, S, MeshU, MeshP, BC, Control);
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
    [fu,fp,fn] = ComputeLoads(BC, MeshU, MeshP, MeshN, Control, QuadU, QuadP);
    
    % analytical solution for 1D case
    if Control.plotansol
        [p_an, u_an] = feval(Control.ansol_type, Control, Material, MeshU, MeshP, BC);
       
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
    [Solution] = SolverTr_UP(Kuu, Kup, Kpp, Kpu, S, fu, fp, BC, Control, Iteration);
    
    % solution in the frequency domain
    if Control.freqDomain
        [SolutionFreq] = SolverFreqTr_UP(phi_u, phi_p, Kuu, Kup, Kpp, S, fu, fp, BC, Control, Iteration);
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
        
        % synthetics
        if Control.fixedDepthPlotON
            Plot.u_synthetic(Control.step+1,:) = Solution.u(Control.ploturow);
            Plot.udot_synthetic(Control.step+1,:) = Solution.udot(Control.ploturow);
            Plot.p_synthetic(Control.step+1,:) = Solution.p(Control.plotprow);
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
    Iteration.fu_old = fu; % load vector
    Iteration.fp_old = fp; % flux vector

    % update variables - frequency domain
    if Control.freqDomain
        Iteration.xuF_old = SolutionFreq.xuF;
        Iteration.xuFdot_old = SolutionFreq.xuFdot;
        Iteration.xpF_old = SolutionFreq.xpF;

        Iteration.uF_old = SolutionFreq.uF; % solid displacement
        Iteration.uFdot_old = SolutionFreq.uFdot; % solid velocity
        Iteration.pF_old = SolutionFreq.pF; % fluid pressure
    end

    % update time and step
    Control.step = Control.step + 1;
end

% plot variables in x for fixed y (2D case)
if Control.fixedDepthPlotON
    Plot.urow = Solution.u(Control.ploturow);
    Plot.udotrow = Solution.udot(Control.ploturow);
    Plot.prow = Solution.p(Control.plotprow);
end