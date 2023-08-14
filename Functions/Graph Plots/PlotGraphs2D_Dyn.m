function PlotGraphs2D_Dyn(MeshU, MeshP, MeshN, Control, Plot, Material, saveGraphs_on)

% initialize figure
figure;
tiledlayout(2,4);

%% pressure vs time
nexttile
plot(Plot.time, Plot.p_time*10^9,'k','LineWidth',2);
hold on
grid on
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

%% solid displacement vs time
nexttile
plot(Plot.time, Plot.u_time,'b','LineWidth',2);
hold on
grid on
xlabel('Time [s]');
ylabel('u (solid) [m]');
title(sprintf('Solid displacement at x = %.2f m, y = %.2f m', MeshU.coords(round(Control.plotu/2),1), MeshU.coords(round(Control.plotu/2),2)));
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
title(sprintf('Solid velocity at x = %.2f m, y = %.2f m', MeshU.coords(round(Control.plotu/2),1), MeshU.coords(round(Control.plotu/2),2)));
% frequency domain solution
if Control.freqDomain
    plot(Plot.time, Plot.uFdot,'m--','LineWidth',2);
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
title(sprintf('Solid acceleration at x = %.2f m, y = %.2f m', MeshU.coords(round(Control.plotu/2),1), MeshU.coords(round(Control.plotu/2),2)));
hold off
if saveGraphs_on
    exportgraphics(gcf,'AccSolid_time.png','Resolution',300)
end

%% porosity vs time
if contains(Control.PMmodel, 'UPN')
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
    %% fluid displacement vs time
    nexttile
    plot(Plot.time, Plot.uf_time,'m','LineWidth',2);
    hold on
    grid on
    xlabel('Time [s]');
    ylabel('u (fluid) [m]');
    title(sprintf('Fluid displacement at x = %.2f m, y = %.2f m', MeshU.coords(round(Control.plotu/2),1), MeshU.coords(round(Control.plotu/2),2)));
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
    title(sprintf('Fluid velocity at x = %.2f m, y = %.2f m', MeshU.coords(round(Control.plotu/2),1), MeshU.coords(round(Control.plotu/2),2)));
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
    title(sprintf('Fluid acceleration at x = %.2f m, y = %.2f m', MeshU.coords(round(Control.plotu/2),1), MeshU.coords(round(Control.plotu/2),2)));
    hold off
    if saveGraphs_on
        exportgraphics(gcf,'AccFluid_time.png','Resolution',300)
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
    grid on
    xlabel('Coordinate [m]');
    ylabel('u [m]');
    title(sprintf('Solid displ. in x at fixed %.0f, %.2f m, t = %.1d s', Control.depthDir, Control.depthplot, Control.tend));
    
    %% velocity in x for fixed coord
    nexttile
    plot(MeshU.coords((Control.ploturow(1:half)+1)./2, Control.depthDir), Plot.udotrow(1:half),'r','LineWidth',2);
    hold on
    grid on
    xlabel('Coordinate [m]');
    ylabel('udot [m]');
    title(sprintf('Solid vel. in x at fixed %.0f, %.2f m, t = %.1d s', Control.depthDir, Control.depthplot, Control.tend));
    
    %% pressure for fixed coord
    nexttile
    plot(MeshP.coords(Control.plotprow, Control.depthDir), Plot.prow*10^9,'k','LineWidth',2);
    hold on
    grid on
    xlabel('Coordinate [m]');
    ylabel('p [Pa]');
    title(sprintf('Press. at fixed %.0f, %.2f m, t = %.1d s', Control.depthDir, Control.depthplot, Control.tend));
    
    %% displacement in y for fixed coord
    nexttile
    plot(MeshU.coords((Control.ploturow(half+1:end))./2, Control.depthDir), Plot.urow(half+1:end),'b','LineWidth',2);
    hold on
    grid on
    xlabel('Coordinate [m]');
    ylabel('u [m]');
    title(sprintf('Solid displ. in y at fixed %.0f, %.2f m, t = %.1d s', Control.depthDir, Control.depthplot, Control.tend));
    
    %% velocity in y for fixed coord
    nexttile
    plot(MeshU.coords((Control.ploturow(half+1:end))./2, Control.depthDir), Plot.udotrow(half+1:end),'r','LineWidth',2);
    hold on
    grid on
    xlabel('Coordinate [m]');
    ylabel('udot [m]');
    title(sprintf('Solid vel. in y at fixed %.0f, %.2f m, t = %.1d s', Control.depthDir, Control.depthplot, Control.tend));
end

if Control.fixedDepthPlotON
    figure;
    %% solid displacement over depth over time (3D plot)
    for t = 1:numel(Plot.time)
        plot3(MeshU.coords((Control.ploturow(half+1:end))./2, Control.depthDir), ones(half,1)*Plot.time(t), Plot.u_synthetic(t,half+1:end), 'k');
        hold on
    end
    xlabel('x [m]');
    ylabel('Time [s]');
    zlabel('u [m]');
    title('Solid displacement in domain over time');
    hold off
    
    % surface plot
    figure;
    surf(MeshU.coords((Control.ploturow(half+1:end))./2, Control.depthDir), Plot.time, Plot.u_synthetic(:,half+1:end));
    hold on
    xlabel('x [m]');
    ylabel('Time [s]');
    zlabel('u [m]');
    title('Solid displacement in domain over time');
    hold off
    
    % waterfall plot
    figure;
    waterfall(MeshU.coords((Control.ploturow(half+1:end))./2, Control.depthDir), Plot.time, Plot.u_synthetic(:,half+1:end));
    hold on
    xlabel('x [m]');
    ylabel('Time [s]');
    zlabel('u [m]');
    title('Solid displacement in domain over time');
    hold off
    
    % waterfall plot
    figure;
    waterfall(MeshP.coords(Control.plotprow, Control.depthDir), Plot.time, Plot.p_synthetic);
    hold on
    xlabel('x [m]');
    ylabel('Time [s]');
    zlabel('p [m]');
    title('Fluid pressure in domain over time');
    hold off

end

end