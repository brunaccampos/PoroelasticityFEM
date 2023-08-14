function PlotGraphs1D_Dyn(Solution, SolutionFreq, MeshU, MeshP, MeshN, Control, Plot, Material, saveGraphs_on)

% initialize figure
figure;
tiledlayout(2,4);

%% pressure vs time
nexttile
plot(Plot.time, Plot.p_time*10^9.*Material.t,'k','LineWidth',2);
hold on
grid on
xlabel('Time [s]');
ylabel('p [Pa]');
title(sprintf('Pressure at x = %.2f m', MeshP.coords(Control.plotp,1)));
% frequency domain solution
if Control.freqDomain
    plot(Plot.time, Plot.pF*10^9.*Material.t,'r--','LineWidth',2);
    legend('Time', 'Frequency');
end
hold off
if saveGraphs_on
    exportgraphics(gcf,'Press_time.png','Resolution',300)
end

%% solid displacement vs time
nexttile
plot(Plot.time, Plot.u_time,'b','LineWidth',2);
hold on
grid on
xlabel('Time [s]');
ylabel('u (solid) [m]');
title(sprintf('Solid displacement at x = %.2f m', MeshU.coords(Control.plotu,1)));
% frequency domain solution
if Control.freqDomain
    plot(Plot.time, Plot.uF,'m--','LineWidth',2);
    legend('Time', 'Frequency');
end
hold off
if saveGraphs_on
    exportgraphics(gcf,'DisplSolid_time.png','Resolution',300)
end

%% solid velocity vs time
nexttile
plot(Plot.time, Plot.udot_time,'r','LineWidth',2);
hold on
grid on
xlabel('Time [s]');
ylabel('udot (solid) [m/s]');
title(sprintf('Solid velocity at x = %.2f m', MeshU.coords(Control.plotu,1)));
% frequency domain solution
if Control.freqDomain
    plot(Plot.time, Plot.uFdot,'g--','LineWidth',2);
    legend('Time', 'Frequency');
end
hold off
if saveGraphs_on
    exportgraphics(gcf,'VelSolid_time.png','Resolution',300)
end

%% solid acceleration vs time
nexttile
plot(Plot.time, Plot.u2dot_time,'Color', [0.9290 0.6940 0.1250], 'LineWidth',2);
hold on
grid on
xlabel('Time [s]');
ylabel('u2dot (solid) [m/s2]');
title(sprintf('Solid acceleration at x = %.2f m', MeshU.coords(Control.plotu,1)));
hold off
if saveGraphs_on
    exportgraphics(gcf,'AccSolid_time.png','Resolution',300)
end

%% pressure vs depth
nexttile
plot(MeshP.coords, Solution.p*10^9.*Material.t,'k','LineWidth',2);
hold on
grid on
xlabel('Column depth [m]');
ylabel('p [Pa]');
title(sprintf('Pressure at t = %.1d s', Control.tend));
% frequency domain solution
if Control.freqDomain
    plot(MeshP.coords, SolutionFreq.pF*10^9.*Material.t,'r--','LineWidth',2);
    legend('Time', 'Frequency');
end
hold off
if saveGraphs_on
    exportgraphics(gcf,'Press_depth.png','Resolution',300)
end

%% solid displacement vs depth
nexttile
plot(MeshU.coords, Solution.u,'b','LineWidth',2);
hold on
grid on
xlabel('Column depth [m]');
ylabel('u (solid) [m]');
title(sprintf('Solid displacement at t = %.1d s', Control.tend));
% frequency domain solution
if Control.freqDomain
    plot(MeshU.coords, SolutionFreq.uF,'m--','LineWidth',2);
    legend('Time', 'Frequency');
end
hold off
if saveGraphs_on
    exportgraphics(gcf,'DisplSolid_depth.png','Resolution',300)
end

%% solid velocity vs depth
nexttile
plot(MeshU.coords, Solution.udot,'r','LineWidth',2);
hold on
grid on
xlabel('Column depth [m]');
ylabel('udot (solid) [m]');
title(sprintf('Solid velocity at t = %.1d s', Control.tend));
% frequency domain solution
if Control.freqDomain
    plot(MeshU.coords, SolutionFreq.uFdot,'g--','LineWidth',2);
    legend('Time', 'Frequency');
end
hold off
if saveGraphs_on
    exportgraphics(gcf,'VelSolid_depth.png','Resolution',300)
end

%% solid acceleration vs depth
nexttile
plot(MeshU.coords, Solution.u2dot,'Color', [0.9290 0.6940 0.1250], 'LineWidth',2);
hold on
grid on
xlabel('Column depth [m]');
ylabel('u2dot (solid) [m]');
title(sprintf('Solid acceleration at t = %.1d s', Control.tend));
hold off
if saveGraphs_on
    exportgraphics(gcf,'AccSolid_depth.png','Resolution',300)
end

if contains(Control.PMmodel, 'UPN')
    figure;
    tiledlayout(1,2);
    %% porosity vs depth
    nexttile
    plot(MeshN.coords, Solution.n ./ Material.n,'g','LineWidth',2);
    hold on
    grid on
    xlabel('Column depth [m]');
    ylabel('Porosity normalized [-]');
    title(sprintf('Porosity norm at t = %.1d s', Control.tend));
    hold off
    if saveGraphs_on
        exportgraphics(gcf,'Poros_depth.png','Resolution',300)
    end

    %% porosity vs time
    nexttile
    plot(Plot.time, Plot.n_time ./ Material.n,'g','LineWidth',2);
    hold on
    grid on
    xlabel('Time [s]');
    ylabel('Porosity normalized [-]');
    title(sprintf('Porosity norm at x = %.2f m', MeshN.coords(Control.plotp,1)));
    hold off
    if saveGraphs_on
        exportgraphics(gcf,'Poros_time.png','Resolution',300)
    end
end

if contains(Control.PMmodel, 'UPU')
    figure;
    tiledlayout(2,3);
    %% fluid displacement vs time
    nexttile
    plot(Plot.time, Plot.uf_time,'m','LineWidth',2);
    hold on
    grid on
    xlabel('Time [s]');
    ylabel('u (fluid) [m]');
    title(sprintf('Fluid displacement at x = %.2f m', MeshU.coords(Control.plotu,1)));
    hold off
    if saveGraphs_on
        exportgraphics(gcf,'DisplFluid_time.png','Resolution',300)
    end
    
    %% fluid velocity vs time
    nexttile
    plot(Plot.time, Plot.ufdot_time,'c','LineWidth',2);
    hold on
    grid on
    xlabel('Time [s]');
    ylabel('udot (fluid) [m/s]');
    title(sprintf('Fluid velocity at x = %.2f m', MeshU.coords(Control.plotu,1)));
    hold off
    if saveGraphs_on
        exportgraphics(gcf,'VelFluid_time.png','Resolution',300)
    end
    
    %% fluid acceleration vs time
    nexttile
    plot(Plot.time, Plot.uf2dot_time,'Color', [0.4660 0.6740 0.1880],'LineWidth',2);
    hold on
    grid on
    xlabel('Time [s]');
    ylabel('u2dot (fluid) [m/s]');
    title(sprintf('Fluid acceleration at x = %.2f m', MeshU.coords(Control.plotu,1)));
    hold off
    if saveGraphs_on
        exportgraphics(gcf,'AccFluid_time.png','Resolution',300)
    end
    
    %% fluid displacement vs depth
    nexttile
    plot(MeshU.coords, Solution.uf,'m','LineWidth',2);
    hold on
    grid on
    xlabel('Column depth [m]');
    ylabel('u (fluid) [m]');
    title(sprintf('Fluid displacement at t = %.1d s', Control.tend));
    hold off
    if saveGraphs_on
        exportgraphics(gcf,'DisplFluid_depth.png','Resolution',300)
    end
    
    %% fluid velocity vs depth
    nexttile
    plot(MeshU.coords, Solution.ufdot,'c','LineWidth',2);
    hold on
    grid on
    xlabel('Column depth [m]');
    ylabel('udot (fluid) [m]');
    title(sprintf('Fluid velocity at t = %.1d s', Control.tend));
    hold off
    if saveGraphs_on
        exportgraphics(gcf,'VelFluid_depth.png','Resolution',300)
    end
    
    %% fluid acceleration vs depth
    nexttile
    plot(MeshU.coords, Solution.uf2dot,'Color', [0.4660 0.6740 0.1880],'LineWidth',2);
    hold on
    grid on
    xlabel('Column depth [m]');
    ylabel('u2dot (fluid) [m]');
    title(sprintf('Fluid acceleration at t = %.1d s', Control.tend));
    hold off
    if saveGraphs_on
        exportgraphics(gcf,'AccFluid_depth.png','Resolution',300)
    end
    
end

if Control.fixedDepthPlotON
    figure;
    %% solid displacement over depth over time (3D plot)
    for t = 1:numel(Plot.time)
        plot3(MeshU.coords, ones(length(MeshU.coords),1)*Plot.time(t), Plot.u_synthetic(t,:), 'k');
        hold on
    end
    xlabel('x [m]');
    ylabel('Time [s]');
    zlabel('u [m]');
    title('Solid displacement in domain over time');
    hold off
    
    % surface plot
    figure;
    surf(MeshU.coords, Plot.time, Plot.u_synthetic);
    hold on
    xlabel('x [m]');
    ylabel('Time [s]');
    zlabel('u [m]');
    title('Solid displacement in domain over time');
    hold off
    
    % waterfall plot
    figure;
    waterfall(MeshU.coords, Plot.time, Plot.u_synthetic);
    hold on
    xlabel('x [m]');
    ylabel('Time [s]');
    zlabel('u [m]');
    title('Solid displacement in domain over time');
    hold off
    
    % waterfall plot
    figure;
    waterfall(MeshP.coords, Plot.time, Plot.p_synthetic);
    hold on
    xlabel('x [m]');
    ylabel('Time [s]');
    zlabel('p [m]');
    title('Fluid pressure in domain over time');
    hold off

end
