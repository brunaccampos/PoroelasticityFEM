function PlotGraphs2D_Tr(MeshU, MeshP, MeshN, Control, Plot, Material, saveGraphs_on)

% initialize figure
figure;
tiledlayout(2,2);

%% pressure vs time
nexttile
plot(Plot.time, Plot.p_time*10^9,'k','LineWidth',2);
hold on
xlabel('Time [s]');
ylabel('p [Pa]');
title(sprintf('Pressure at x = %.2f m, y = %.2f m', MeshP.coords(round(Control.plotp/2),1), MeshP.coords(round(Control.plotp/2),2)));
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
    exportgraphics(gcf,'Displ_time.png','Resolution',300)
end
%% solid velocity vs time
nexttile
plot(Plot.time, Plot.udot_time,'b','LineWidth',2);
hold on
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
    exportgraphics(gcf,'Vel_time.png','Resolution',300)
end
%% porosity vs time
if contains(Control.PMmodel, 'UPN')
    nexttile
    plot(Plot.time, Plot.n_time ./ Material.n,'g','LineWidth',2);
    hold on
    xlabel('Time [s]');
    ylabel('Porosity normalized [-]');
    title(sprintf('Porosity norm at x = %.2f m', MeshN.coords(round(Control.plotp/2),1), MeshN.coords(round(Control.plotp/2),2)));
    hold off
    if saveGraphs_on
        exportgraphics(gcf,'Poros_time.png','Resolution',300)
    end
end

% initialize figure
figure;
tiledlayout(2,3);

if isfield(Control, 'depthplot')
    %% displacement for fixed coord
    nexttile
    plot(MeshU.coords(Control.ploturow./Control.DOFplot, Control.depthDir), Plot.urow,'b','LineWidth',2);
    hold on
    xlabel('Coordinate [m]');
    ylabel('u [m]');
    title(sprintf('Solid displ. at fixed %.2f m, DOF %.0f, t = %.0f s', Control.depthplot, Control.DOFplot, Control.tend));
    %% velocity for fixed coord
    nexttile
    plot(MeshU.coords(Control.ploturow./Control.DOFplot, Control.depthDir), Plot.udotrow,'b','LineWidth',2);
    hold on
    xlabel('Coordinate [m]');
    ylabel('udot [m]');
    title(sprintf('Solid vel. at fixed %.2f m, DOF %.0f, t = %.0f s', Control.depthplot, Control.DOFplot, Control.tend));
    %% pressure for fixed coord
    nexttile
    plot(MeshP.coords(Control.plotprow,2), Plot.prow*10^9,'b','LineWidth',2);
    hold on
    xlabel('Coordinate [m]');
    ylabel('p [Pa]');
    title(sprintf('Press. at fixed %.2f m, t = %.0f s', Control.depthplot, Control.tend));
    
    if contains (Control.PMmodel,'UPN')
        %% porosity for fixed coord
        nexttile
        plot(MeshN.coords(Control.plotprow,2), (Plot.nrow - Material.n)./ Material.n,'g','LineWidth',2);
        hold on
        xlabel('y [m]');
        ylabel('Change in porosity normalized [-]');
        title(sprintf('Porosity at fixed %.2f m, t = %.0f s', Control.depthplot, Control.tend));
    end
end

end