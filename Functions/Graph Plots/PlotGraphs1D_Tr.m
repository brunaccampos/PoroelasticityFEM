function PlotGraphs1D_Tr(Solution, SolutionFreq, Material, MeshU, MeshP, MeshN, Control, Plot, saveGraphs_on)

% initialize figure
figure;
tiledlayout(2,3);

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
% analytical solution
if Control.plotansol
    plot(Plot.time, Plot.pan_time*10^9.*Material.t,'k:','LineWidth',2);
    legend('Numerical', 'Analytical');
end
% adapt legend if both frequency and analytical solutions are plotted
if Control.freqDomain && Control.plotansol
    legend('Time', 'Frequency', 'Analytical');
end
hold off
if saveGraphs_on
    exportgraphics(gcf,'Press_time.png','Resolution',300)
end

%% solid displacement vs time
nexttile
plot(Plot.time, Plot.u_time, 'b', 'LineWidth',2);
hold on
grid on
xlabel('Time [s]');
ylabel('u (solid) [m]');
title(sprintf('Solid displacement at x = %.2f m', MeshU.coords(Control.plotu,1)));
% frequency domain solution
if Control.freqDomain
    plot(Plot.time, Plot.uF, 'm--', 'LineWidth', 2);
    legend('Time', 'Frequency');
end
% analytical solution
if Control.plotansol
    plot(Plot.time, Plot.uan_time,'k:','LineWidth',2);
    legend('Numerical', 'Analytical');
end
% adapt legend if both frequency and analytical solutions are plotted
if Control.freqDomain && Control.plotansol
    legend('Time', 'Frequency', 'Analytical');
end
hold off
if saveGraphs_on
    exportgraphics(gcf,'Displ_time.png','Resolution',300)
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
    plot(Plot.time, Plot.uFdot, 'g--', 'LineWidth',2);
    legend('Time', 'Frequency');
end
hold off
if saveGraphs_on
    exportgraphics(gcf,'Vel_time.png','Resolution',300)
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
% analytical solution
if Control.plotansol
    plot(MeshP.coords, Plot.pan_space*10^9.*Material.t,'k:','LineWidth',2);
    legend('Numerical', 'Analytical');
end
% adapt legend if both frequency and analytical solutions are plotted
if Control.freqDomain && Control.plotansol
    legend('Time', 'Frequency', 'Analytical');
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
% analytical solution
if Control.plotansol
    plot(MeshU.coords, Plot.uan_space,'k:','LineWidth',2);
    legend('Numerical', 'Analytical');
end
% adapt legend if both frequency and analytical solutions are plotted
if Control.freqDomain && Control.plotansol
    legend('Time', 'Frequency', 'Analytical');
end
hold off
if saveGraphs_on
    exportgraphics(gcf,'Displ_depth.png','Resolution',300)
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
    exportgraphics(gcf,'Displ_depth.png','Resolution',300)
end

if contains(Control.PMmodel, 'UPN')
    figure;
    tiledlayout(1,2);
    %% porosity vs depth
    nexttile
    semilogx(MeshN.coords, (Solution.n-Material.eta0) ./ Material.eta0,'g','LineWidth',2);
    grid on
    hold on
    xlabel('Column depth [m]');
    ylabel('Change in porosity normalized [-]');
    title(sprintf('Porosity norm at t = %.1d s', Control.tend));
    hold off
    if saveGraphs_on
        exportgraphics(gcf,'Poros_depth.png','Resolution',300)
    end
    
    %% porosity vs time
    nexttile
    semilogx(Plot.time, (Plot.n_time-Material.eta0) ./ Material.eta0,'g','LineWidth',2);
    hold on
    grid on
    xlabel('Time [s]');
    ylabel('Change in porosity normalized [-]');
    title(sprintf('Porosity norm at x = %.2f m', MeshN.coords(Control.plotp,1)));
    hold off
    if saveGraphs_on
        exportgraphics(gcf,'Poros_time.png','Resolution',300)
    end
end

if contains(Control.PMmodel, 'UPU') || contains(Control.PMmodel, 'UPV') || contains(Control.PMmodel, 'UPW')
    figure;
    tiledlayout(1,2);
    %% porosity vs depth
    nexttile;
    semilogx(MeshU.coords, (Solution.n-Material.eta0) ./ Material.eta0,'g','LineWidth',2);
    grid on
    hold on
    xlabel('Column depth [m]');
    ylabel('Change in porosity normalized [-]');
    title(sprintf('Porosity norm at t = %.1d s', Control.tend));
    hold off
    if saveGraphs_on
        exportgraphics(gcf,'Poros_depth.png','Resolution',300)
    end
    %% porosity gradient vs depth
    nexttile;
    plot(MeshU.coords, Solution.ndot,'g','LineWidth',2);
    grid on
    hold on
    xlabel('Column depth [m]');
    ylabel('Porosity gradient [-]');
    title(sprintf('Porosity gradient (dt) at t = %.1d s', Control.tend));
    hold off
    if saveGraphs_on
        exportgraphics(gcf,'PorosDot_depth.png','Resolution',300)
    end
end

if Control.fixedDepthPlotON
    %% solid displacement over depth over time (3D plot)
    figure;
    waterfall(MeshU.coords, Plot.time, Plot.u_synthetic.*10^7);
    hold on
    xlabel('x (m)','interpreter','latex');
    ylabel('Time (s)','interpreter','latex');
    zlabel('Displacement ($\times 10^{-7}$ m)','interpreter','latex');
    title('Solid displacement in domain over time','interpreter','latex');
    view(135, 30);
    hold off
    exportgraphics(gcf,'SolidDisp_3D.png','Resolution',300);
    savefig(gcf,'SolidDisp_3D.fig');

    %% fluid pressure over depth over time (3D plot)
    figure;
    waterfall(MeshP.coords, Plot.time, Plot.p_synthetic.*10^9);
    hold on
    xlabel('x (m)','interpreter','latex');
    ylabel('Time (s)','interpreter','latex');
    zlabel('Pressure (Pa)','interpreter','latex');
    title('Fluid pressure in domain over time','interpreter','latex');
    view(135, 30);
    hold off
    exportgraphics(gcf,'Pressure_3D.png','Resolution',300);
    savefig(gcf,'Pressure_3D.fig');
end

end