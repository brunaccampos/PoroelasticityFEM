function PlotGraphs(Solution, SolutionFreq, Material, MeshU, MeshP, MeshN, Control, Plot, saveGraphs_on)
% ------------------------------------------------------------------------
% Plot results in graphs
% ------------------------------------------------------------------------
% Input: Mesh, Control, Plot, u, p
% ------------------------------------------------------------------------
% Output: plots of
%               displacement vs depth
%               pressure vs depth
%               displacement vs time at specific node
%               pressure vs time at specific node
% ------------------------------------------------------------------------

if Control.steady
    if MeshU.nsd == 1
        plot1Dsteady(Solution, SolutionFreq, Material, MeshU, MeshP, MeshN, Control, Plot, saveGraphs_on);
    else
        plot2Dsteady(Solution, MeshU, MeshP, Control, Plot, saveGraphs_on);
    end
else
    if MeshU.nsd == 1
        plot1Ddynamic(Solution, SolutionFreq, MeshU, MeshP, MeshN, Control, Plot, Material, saveGraphs_on);
    else
        plot2Ddynamic(Solution, MeshU, MeshP, Control, Plot, saveGraphs_on);
    end
end

end

% ------------------------------------------------------------------------
% PLOTS FOR 1D QUASI-STEADY CASE
% ------------------------------------------------------------------------
function plot1Dsteady(Solution, SolutionFreq, Material, MeshU, MeshP, MeshN, Control, Plot, saveGraphs_on)
%% pressure vs depth
figure;
plot(MeshP.coords, Solution.p*10^9,'k','LineWidth',2);
hold on
xlabel('Column depth [m]');
ylabel('p [Pa]');
title(sprintf('Pressure at t = %.0f s', Control.tend));
% frequency domain solution
if Control.freqDomain
    plot(MeshP.coords, SolutionFreq.pF*10^9,'r--','LineWidth',2);
    legend('Time', 'Frequency');
end
% analytical solution
if Control.plotansol
    plot(MeshP.coords, Plot.pan_space*10^9,'k:','LineWidth',2);
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
figure;
plot(Plot.time, Plot.p_time*10^9,'k','LineWidth',2);
hold on
xlabel('Time [s]');
ylabel('p [Pa]');
title(sprintf('Pressure at x = %.2f m', MeshP.coords(Control.plotp,1)));
% frequency domain solution
if Control.freqDomain
    plot(Plot.time, Plot.pF*10^9,'r--','LineWidth',2);
    legend('Time', 'Frequency');
end
% analytical solution
if Control.plotansol
    plot(Plot.time, Plot.pan_time*10^9,'k:','LineWidth',2);
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
%% displacement vs depth
figure;
plot(MeshU.coords, Solution.u,'b','LineWidth',2);
hold on
xlabel('Column depth [m]');
ylabel('u [m]');
title(sprintf('Displacement at t = %.0f s', Control.tend));
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
%% displacement vs time
figure;
plot(Plot.time, Plot.u_time, 'b', 'LineWidth',2);
hold on
xlabel('Time [s]');
ylabel('u [m]');
title(sprintf('Solid skeleton displacement at x = %.2f m', MeshU.coords(Control.plotu,1)));
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
%% velocity vs time
figure;
plot(Plot.time, Plot.udot_time,'b','LineWidth',2);
hold on
xlabel('Time [s]');
ylabel('udot [m/s]');
title(sprintf('Solid skeleton velocity at x = %.2f m', MeshU.coords(Control.plotu,1)));
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
if ~Control.Biotmodel
    figure;
    plot(MeshN.coords, Solution.n ./ Material.n,'g','LineWidth',2);
    hold on
    xlabel('Column depth [m]');
    ylabel('Porosity normalized [-]');
    title(sprintf('Porosity norm at t = %.0f s', Control.tend));
    hold off
    if saveGraphs_on
        exportgraphics(gcf,'Poros_depth.png','Resolution',300)
    end
end
%% porosity vs time
if ~Control.Biotmodel
    figure;
    plot(Plot.time, Plot.n_time ./ Material.n,'g','LineWidth',2);
    hold on
    xlabel('Time [s]');
    ylabel('Porosity normalized [-]');
    title(sprintf('Porosity norm at x = %.2f m', MeshN.coords(Control.plotp,1)));
    hold off
    if saveGraphs_on
        exportgraphics(gcf,'Poros_time.png','Resolution',300)
    end
end
end

% ------------------------------------------------------------------------
% PLOTS FOR 2D QUASI-STEADY CASE
% ------------------------------------------------------------------------
function plot2Dsteady(Solution, MeshU, MeshP, Control, Plot, saveGraphs_on)
%% pressure vs time
figure;
plot(Plot.time, Plot.p_time*10^9,'k','LineWidth',2);
hold on
xlabel('Time [s]');
ylabel('p [Pa]');
title(sprintf('Pressure at x = %.2f m, y = %.2f m', MeshP.coords(Control.plotp,1), MeshP.coords(Control.plotp,2)));
% frequency domain solution
if Control.freqDomain
    plot(Plot.time, Plot.pF*10^9,'r--','LineWidth',2);
    legend('Time', 'Frequency');
end
hold off
if saveGraphs_on
    exportgraphics(gcf,'Press_time.png','Resolution',300)
end
%% displacement vs time
figure;
plot(Plot.time, Plot.u_time,'b','LineWidth',2);
hold on
xlabel('Time [s]');
ylabel('u [m]');
title(sprintf('Solid skeleton displacement at x = %.2f m, y = %.2f m', MeshU.coords(round(Control.plotu/2),1), MeshU.coords(round(Control.plotu/2),2)));
% frequency domain solution
if Control.freqDomain
    plot(Plot.time, Plot.uF,'m--','LineWidth',2);
    legend('Time', 'Frequency');
end
hold off
if saveGraphs_on
    exportgraphics(gcf,'Displ_time.png','Resolution',300)
end
%% velocity vs time
figure;
plot(Plot.time, Plot.udot_time,'b','LineWidth',2);
hold on
xlabel('Time [s]');
ylabel('udot [m/s]');
title(sprintf('Solid skeleton velocity at x = %.2f m, y = %.2f m', MeshU.coords(round(Control.plotu/2),1), MeshU.coords(round(Control.plotu/2),2)));
% frequency domain solution
if Control.freqDomain
    plot(Plot.time, Plot.uFdot,'m--','LineWidth',2);
    legend('Time', 'Frequency');
end
hold off
if saveGraphs_on
    exportgraphics(gcf,'Vel_time.png','Resolution',300)
end
% %% pressure over mesh
% figure;
% trisurf(MeshP.conn, MeshP.coords(:,1), MeshP.coords(:,2), Solution.p*1e9);
% view(2)
% hold on
% colorbar
% xlabel('x [m]');
% ylabel('y [m]');
% title(sprintf('Pressure [Pa] at t = %.0f s', Control.tend));
% hold off
end

% ------------------------------------------------------------------------
% PLOTS FOR 1D DYNAMIC CASE
% ------------------------------------------------------------------------
function plot1Ddynamic(Solution, SolutionFreq, MeshU, MeshP, MeshN, Control, Plot, Material, saveGraphs_on)
%% pressure vs depth
figure;
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
hold off
if saveGraphs_on
    exportgraphics(gcf,'Press_depth.png','Resolution',300)
end
%% pressure vs time
figure;
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
hold off
if saveGraphs_on
    exportgraphics(gcf,'Press_time.png','Resolution',300)
end
%% displacement vs depth
figure;
plot(MeshU.coords, Solution.u,'b','LineWidth',2);
hold on
xlabel('Column depth [m]');
ylabel('u [m]');
title(sprintf('Displacement at t = %.0f s', Control.tend));
% frequency domain solution
if Control.freqDomain
    plot(MeshU.coords, SolutionFreq.uF,'m--','LineWidth',2);
    legend('Time', 'Frequency');
end
hold off
if saveGraphs_on
    exportgraphics(gcf,'Displ_depth.png','Resolution',300)
end
%% displacement vs time
figure;
plot(Plot.time, Plot.u_time,'b','LineWidth',2);
hold on
xlabel('Time [s]');
ylabel('u [m]');
title(sprintf('Solid skeleton displacement at x = %.2f m', MeshU.coords(Control.plotu,1)));
% frequency domain solution
if Control.freqDomain
    plot(Plot.time, Plot.uF,'m--','LineWidth',2);
    legend('Time', 'Frequency');
end
hold off
if saveGraphs_on
    exportgraphics(gcf,'Displ_time.png','Resolution',300)
end
%% velocity vs time
figure;
plot(Plot.time, Plot.udot_time,'b','LineWidth',2);
hold on
xlabel('Time [s]');
ylabel('udot [m/s]');
title(sprintf('Solid skeleton velocity at x = %.2f m', MeshU.coords(Control.plotu,1)));
% frequency domain solution
if Control.freqDomain
    plot(Plot.time, Plot.uFdot,'m--','LineWidth',2);
    legend('Time', 'Frequency');
end
hold off
if saveGraphs_on
    exportgraphics(gcf,'Vel_time.png','Resolution',300)
end
%% porosity vs depth
if ~Control.Biotmodel
    figure;
    plot(MeshN.coords, Solution.n ./ Material.n,'g','LineWidth',2);
    hold on
    xlabel('Column depth [m]');
    ylabel('Porosity normalized [-]');
    title(sprintf('Porosity norm at t = %.0f s', Control.tend));
    hold off
    if saveGraphs_on
        exportgraphics(gcf,'Poros_depth.png','Resolution',300)
    end
end
%% porosity vs time
if ~Control.Biotmodel
    figure;
    plot(Plot.time, Plot.n_time ./ Material.n,'g','LineWidth',2);
    hold on
    xlabel('Time [s]');
    ylabel('Porosity normalized [-]');
    title(sprintf('Porosity norm at x = %.2f m', MeshN.coords(Control.plotp,1)));
    hold off
    if saveGraphs_on
        exportgraphics(gcf,'Poros_time.png','Resolution',300)
    end
end

end

% ------------------------------------------------------------------------
% PLOTS FOR 2D DYNAMIC CASE
% ------------------------------------------------------------------------
function plot2Ddynamic(Solution, MeshU, MeshP, Control, Plot, saveGraphs_on)
%% pressure vs time
figure;
plot(Plot.time, Plot.p_time*10^9,'k','LineWidth',2);
hold on
xlabel('Time [s]');
ylabel('p [Pa]');
title(sprintf('Pressure at x = %.2f m, y = %.2f m', MeshP.coords(Control.plotp,1), MeshP.coords(Control.plotp,2)));
% frequency domain solution
if Control.freqDomain
    plot(Plot.time, Plot.pF*10^9,'r--','LineWidth',2);
    legend('Time', 'Frequency');
end
hold off
if saveGraphs_on
    exportgraphics(gcf,'Press_time.png','Resolution',300)
end
%% displacement vs time
figure;
plot(Plot.time, Plot.u_time,'b','LineWidth',2);
hold on
xlabel('Time [s]');
ylabel('u [m]');
title(sprintf('Solid skeleton displacement at x = %.2f m, y = %.2f m', MeshU.coords(round(Control.plotu/2),1), MeshU.coords(round(Control.plotu/2),2)));
% frequency domain solution
if Control.freqDomain
    plot(Plot.time, Plot.uF,'m--','LineWidth',2);
    legend('Time', 'Frequency');
end
hold off
if saveGraphs_on
    exportgraphics(gcf,'Displ_time.png','Resolution',300)
end
%% velocity vs time
figure;
plot(Plot.time, Plot.udot_time,'b','LineWidth',2);
hold on
xlabel('Time [s]');
ylabel('udot [m/s]');
title(sprintf('Solid skeleton velocity at x = %.2f m, y = %.2f m', MeshU.coords(round(Control.plotu/2),1), MeshU.coords(round(Control.plotu/2),2)));
% frequency domain solution
if Control.freqDomain
    plot(Plot.time, Plot.uFdot,'m--','LineWidth',2);
    legend('Time', 'Frequency');
end
hold off
if saveGraphs_on
    exportgraphics(gcf,'Vel_time.png','Resolution',300)
end
% %% pressure over mesh
% figure;
% trisurf(MeshP.conn, MeshP.coords(:,1), MeshP.coords(:,2), Solution.p*1e9);
% view(2)
% hold on
% colorbar
% xlabel('x [m]');
% ylabel('y [m]');
% title(sprintf('Pressure [Pa] at t = %.0f s', Control.tend));
% hold off

end