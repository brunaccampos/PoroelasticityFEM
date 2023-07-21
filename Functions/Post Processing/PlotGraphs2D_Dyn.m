function PlotGraphs2D_Dyn(MeshU, MeshP, MeshN, Control, Plot, Material, saveGraphs_on)

% initialize figure
figure;
tiledlayout(2,2);

%% pressure vs time
nexttile
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
nexttile
plot(Plot.time, Plot.u_time,'b','LineWidth',2);
hold on
xlabel('Time [s]');
ylabel('u [m]');
title(sprintf('Solid displacement at x = %.2f m, y = %.2f m', MeshU.coords(round(Control.plotu/2),1), MeshU.coords(round(Control.plotu/2),2)));
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
nexttile
plot(Plot.time, Plot.udot_time,'r','LineWidth',2);
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

%% porosity vs time
if contains(Control.PMmodel, 'UPN')
    nexttile
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

if Control.fixedDepthPlotON
    % initialize figure
    figure;
    tiledlayout(2,3);

    % find half of array
    half = ceil(length(Control.ploturow)/2);
    
    %% displacement in x for fixed coord
    nexttile
    plot(MeshU.coords((Control.ploturow(1:half)+1)./2, Control.depthDir), Plot.urow(1:half),'b','LineWidth',2);
    hold on
    xlabel('Coordinate [m]');
    ylabel('u [m]');
    title(sprintf('Solid displ. in x at fixed %.0f, %.2f m, t = %.1d s', Control.depthDir, Control.depthplot, Control.tend));
    
    %% velocity in x for fixed coord
    nexttile
    plot(MeshU.coords((Control.ploturow(1:half)+1)./2, Control.depthDir), Plot.udotrow(1:half),'r','LineWidth',2);
    hold on
    xlabel('Coordinate [m]');
    ylabel('udot [m]');
    title(sprintf('Solid vel. in x at fixed %.0f, %.2f m, t = %.1d s', Control.depthDir, Control.depthplot, Control.tend));
    
    %% pressure for fixed coord
    nexttile
    plot(MeshP.coords(Control.plotprow, Control.depthDir), Plot.prow*10^9,'k','LineWidth',2);
    hold on
    xlabel('Coordinate [m]');
    ylabel('p [Pa]');
    title(sprintf('Press. at fixed %.0f, %.2f m, t = %.1d s', Control.depthDir, Control.depthplot, Control.tend));
    
    %% displacement in y for fixed coord
    nexttile
    plot(MeshU.coords((Control.ploturow(half+1:end))./2, Control.depthDir), Plot.urow(half+1:end),'b','LineWidth',2);
    hold on
    xlabel('Coordinate [m]');
    ylabel('u [m]');
    title(sprintf('Solid displ. in y at fixed %.0f, %.2f m, t = %.1d s', Control.depthDir, Control.depthplot, Control.tend));
    
    %% velocity in y for fixed coord
    nexttile
    plot(MeshU.coords((Control.ploturow(half+1:end))./2, Control.depthDir), Plot.udotrow(half+1:end),'r','LineWidth',2);
    hold on
    xlabel('Coordinate [m]');
    ylabel('udot [m]');
    title(sprintf('Solid vel. in y at fixed %.0f, %.2f m, t = %.1d s', Control.depthDir, Control.depthplot, Control.tend));
end

end