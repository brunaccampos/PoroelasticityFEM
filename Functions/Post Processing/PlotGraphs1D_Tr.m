function PlotGraphs1D_Tr(Solution, SolutionFreq, Material, MeshU, MeshP, MeshN, Control, Plot, saveGraphs_on)

% initialize figure
figure;
tiledlayout(2,4);

%% pressure vs depth
nexttile
plot(MeshP.coords, Solution.p*10^9.*Material.t,'k','LineWidth',2);
hold on
xlabel('Column depth [m]');
ylabel('p [Pa]');
title(sprintf('Pressure at t = %.0f s', Control.tend));
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

%% pressure vs time
nexttile
plot(Plot.time, Plot.p_time*10^9.*Material.t,'k','LineWidth',2);
hold on
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

%% solid displacement vs depth
nexttile
plot(MeshU.coords, Solution.u,'b','LineWidth',2);
hold on
xlabel('Column depth [m]');
ylabel('u (solid) [m]');
title(sprintf('Solid displacement at t = %.0f s', Control.tend));
% frequency domain solution
if Control.freqDomain
    plot(MeshU.coords, SolutionFreq.uF,'m--','LineWidth',2);
    legend('Time', 'Frequency');
end
% analytical solution
if Control.plotansol
    %     syms x
    %     fplot((-x^2 +2*x), [0,1], 'LineWidth',2);
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

%% solid displacement vs time
nexttile
plot(Plot.time, Plot.u_time, 'b', 'LineWidth',2);
hold on
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
plot(Plot.time, Plot.udot_time,'b','LineWidth',2);
hold on
xlabel('Time [s]');
ylabel('udot (solid) [m/s]');
title(sprintf('Solid velocity at x = %.2f m', MeshU.coords(Control.plotu,1)));
% frequency domain solution
if Control.freqDomain
    plot(Plot.time, Plot.uFdot, 'm--', 'LineWidth',2);
    legend('Time', 'Frequency');
end
hold off
if saveGraphs_on
    exportgraphics(gcf,'Vel_time.png','Resolution',300)
end

%% porosity vs depth
if contains(Control.PMmodel, 'UPN')
    nexttile
    semilogx(MeshN.coords, (Solution.n-Material.n) ./ Material.n,'g','LineWidth',2);
    hold on
    xlabel('Column depth [m]');
    ylabel('Change in porosity normalized [-]');
    title(sprintf('Porosity norm at t = %.0f s', Control.tend));
    hold off
    if saveGraphs_on
        exportgraphics(gcf,'Poros_depth.png','Resolution',300)
    end
    
%% porosity vs time
    nexttile
    semilogx(Plot.time, (Plot.n_time-Material.n) ./ Material.n,'g','LineWidth',2);
    hold on
    xlabel('Time [s]');
    ylabel('Change in porosity normalized [-]');
    title(sprintf('Porosity norm at x = %.2f m', MeshN.coords(Control.plotp,1)));
    hold off
    if saveGraphs_on
        exportgraphics(gcf,'Poros_time.png','Resolution',300)
    end
end

end