% Porous Media Simulation
% ------------------------------------------------------------------------
% Created by Bruna Campos
% bccampos@uwaterloo.ca
% Department of Civil Engineering, University of Waterloo
% January 2022
% ------------------------------------------------------------------------
% Reference: https://github.com/GCMLab
% ------------------------------------------------------------------------
% version 3 (11/10/22): cleaning unecessary code
% ------------------------------------------------------------------------

%% Delete past files
if plot2vtk || plot2csv_on
    if progress_on
        fprintf('%.2f: Deleting old vtk and csv files...\n', toc);
    end
    sol = fullfile(vtk_dir, '*');
    delete(sol)
end

%% Input problem
if progress_on
    fprintf('%.2f: Reading config file...\n', toc);
end
[Material, MeshU, MeshP, MeshN, BC, Control] = feval(config_name, ConfigDir, progress_on);

%% Quadrature points
Quad = GlobalQuad(MeshU, Control);

%% Solution parameters
Control.t = 0;  % initial simulation time
Control.step = 0; % step counter

%% Assemble system matrices
disp([num2str(toc),': Assembling System Matrices...']);

if Control.steady %% quasi-steady/transient case
    if Control.Biotmodel % Biot
        [Kuu, Kup, Kpp, S] = ComputeSystemMatrices_Transient(Material, MeshU, MeshP, Quad);
    else % Spanos
        [Kuu, Kup, Kpp, Kpu, S, Kpn, Knn, Knu, Knp, Kun] = ComputeSystemMatrices_Transient_v3(Material, MeshU, MeshP, MeshN,Quad);
    end
else %% dynamic case
    if Control.Biotmodel % Biot
        [Kuu, Kup, Kpp, M, Mhat, S] = ComputeSystemMatrices_Dynamic(Material, MeshU, MeshP, Quad);
    else % Spanos
        [Muu, Mpu, Mnu, Kuu, Kup, Kpp, Kpu, S, Kpn, Knn, Knu, Knp, Kun] = ComputeSystemMatrices_Dynamic_v3(Material, MeshU, MeshP, MeshN, Quad);
    end
end

%% Assemble system load vectors
[fu,fp,fn] = ComputeSystemLoads(BC, MeshU, MeshP, MeshN, Control, Quad);

%% Solve eigenproblem
if Control.freqDomain
    disp([num2str(toc),': Solving Uncoupled Eigenproblems...']);   
    if Control.steady %% quasi-steady/transient case
        [phi_u, omega2_u, phi_p, omega2_p] = SolveEig_Transient(Kuu, Kup, Kpp, MeshU, MeshP, BC, Control);
    else %% dynamic case
        [phi_u, omega2_u, phi_p, omega2_p] = SolveEig_Dynamic(Kuu, Kup, Kpp, M, Mhat, MeshU, MeshP, BC, Control);        
    end
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
    
    if ~Control.Biotmodel
        Solution.n = Iteration.n_old;
    else
        Solution.n = [];
    end
    % save vtk file
    PostProcessing(Solution, Material, MeshU, MeshP, MeshN, Control, BC, config_name, vtk_dir);
    % update time step
    Control.step = 1;
end

%% Solve system
while Control.t < Control.tend
    fprintf('\n Step %d \n', Control.step);
    
    %%%% change for this load computation if time dependent
%     [fu,fp,fn] = ComputeSystemLoads(BC, MeshU, MeshP, MeshN, Control, Quad);
    
    % linear solver
    if Control.steady
        if Control.plotansol
            % analytical solution for 1D case and 1/Q = 0 (null S matrix)
            % (incompressible solid and fluid)
%             [~, p_an, u_an] = getAnalyticResult(Material, MeshU, MeshP, BC, Control);
            [p_an, u_an] = getAnalyticResult_v2(Material, MeshU, MeshP, BC, Control);

            % store variables over length
            Plot.pan_space = p_an;
            Plot.uan_space = u_an;
            
            % store variables over time
            if Control.step < length(Plot.time)
                Plot.uan_time(Control.step+1,:) = u_an(Control.plotu, 1);
                Plot.pan_time(Control.step+1,:) = p_an(Control.plotp, 1);
            end
        end
        % quasi-steady/transient case
        if Control.Biotmodel
            [Solution] = SolverTransient_v4(Kuu, Kup, Kpp, S, fu, fp, BC, Control, Iteration);
            if Control.freqDomain
                [SolutionFreq] = SolverTransientFreq(phi_u, omega2_u, phi_p, omega2_p, Kuu, Kup, Kpp, S, fu, fp, BC, Control, Iteration);
            end
        else
            [Solution] = SolverTransient_v5(Kuu, Kup, Kpp, Kpu, S, Kpn, Knn, Knu, Knp, Kun, fu, fp, fn, BC, Control, Iteration);
        end
    else
        % dynamic case
        if Control.Biotmodel
            [Solution] = SolverDynamic_v4(Kuu, Kup, Kpp, M, Mhat, S, fu, fp, BC, Control, Iteration);
            if Control.freqDomain
                [SolutionFreq] = SolverDynamicFreq(phi_u, omega2_u, phi_p, omega2_p, Kuu, Kup, Kpp, M, Mhat, S, fu, fp, BC, Control, Iteration);
            end
        else
            [Solution] = SolverDynamic_v6(Muu, Mpu, Mnu, Kuu, Kup, Kpp, Kpu, S, Kpn, Knn, Knu, Knp, Kun, fu, fp, fn, BC, Control, Iteration);
        end
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
        
        % plot porosity vs time
        if ~Control.Biotmodel
            Plot.n_time(Control.step+1,:) = Solution.n(Control.plotp, 1) ;
        end
        
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
    Iteration.fu_old = fu; % load vector
    Iteration.fp_old = fp; % flux vector
    
    if ~Control.Biotmodel
        Iteration.n_old = Solution.n; % medium porosity
    end
    
    % update variables - frequency domain
    if Control.freqDomain
        Iteration.xuF_old = SolutionFreq.xuF;
        Iteration.xuFdot_old = SolutionFreq.xuFdot;
        Iteration.xpF_old = SolutionFreq.xpF;
        
        Iteration.uF_old = SolutionFreq.uF; % solid displacement
        Iteration.uFdot_old = SolutionFreq.uFdot; % solid velocity
        Iteration.pF_old = SolutionFreq.pF; % fluid pressure
    end
    
    if ~Control.steady
        % time domain
        Iteration.u2dot_old = Solution.u2dot; % solid acceleration
        Iteration.pdot_old = Solution.pdot; % fluid pressure gradient
        if ~Control.Biotmodel
            Iteration.ndot_old = Solution.ndot; % porosity variation
        end
        
        % frequency domain
        if Control.freqDomain
            Iteration.xuF2dot_old = SolutionFreq.xuF2dot;
            Iteration.xpFdot_old = SolutionFreq.xpFdot;
            
            Iteration.uF2dot_old = SolutionFreq.uF2dot; % solid acceleration
            Iteration.pFdot_old = SolutionFreq.pFdot; % fluid pressure gradient
        end
    end
    
    % update time and step
    Control.t = Control.t + Control.dt;
    Control.step = Control.step + 1;
end

%% Post processing
% plot figures
if Control.freqDomain
    PlotGraphs(Solution, SolutionFreq, Material, MeshU, MeshP, MeshN, Control, Plot, saveGraphs_on);
else
    PlotGraphs(Solution, [], Material, MeshU, MeshP, MeshN, Control, Plot, saveGraphs_on);
end

% plot natural frequencies
if Control.freqDomain
    PlotModeShapes(phi_u, omega2_u, phi_p, omega2_p, MeshU, MeshP, Control, BC, config_name, vtk_dir);
end

% export CSV file
if plot2csv_on
    Plot2csv(config_name, vtk_dir, Plot);
end

% store results for error computation
u = Solution.u;
p = Solution.p;
q = ComputeFluidFlux(Material, MeshP, p);
u_an = Plot.uan_space;
p_an = Plot.pan_space;
q_an = ComputeFluidFlux(Material, MeshP, p_an);
save('Results.mat', 'Control', 'Plot', 'Quad', 'MeshU', 'u', 'u_an', 'MeshP', 'p', 'p_an', 'q', 'q_an');
