function PlotGraphs1D_Dyn(Solution, SolutionFreq, MeshU, MeshP, MeshN, Control, Plot, Material, saveGraphs_on)

figure;
%% pressure vs depth
subplot(2,3,1);
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
subplot(2,3,4);
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
%% solid displacement vs depth
subplot(2,3,2);
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
hold off
if saveGraphs_on
    exportgraphics(gcf,'Displ_depth.png','Resolution',300)
end
%% solid displacement vs time
subplot(2,3,5);
plot(Plot.time, Plot.u_time,'b','LineWidth',2);
hold on
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
    exportgraphics(gcf,'Displ_time.png','Resolution',300)
end
%% solid velocity vs time
subplot(2,3,3);
plot(Plot.time, Plot.udot_time,'r','LineWidth',2);
hold on
xlabel('Time [s]');
ylabel('udot (solid) [m/s]');
title(sprintf('Solid velocity at x = %.2f m', MeshU.coords(Control.plotu,1)));
% frequency domain solution
if Control.freqDomain
    plot(Plot.time, Plot.uFdot,'m--','LineWidth',2);
    legend('Time', 'Frequency');
end
hold off
if saveGraphs_on
    exportgraphics(gcf,'Vel_time.png','Resolution',300)
end
%% solid acceleration vs time
subplot(2,3,6);
plot(Plot.time, Plot.u2dot_time,'r','LineWidth',2);
hold on
xlabel('Time [s]');
ylabel('u2dot (solid) [m/s]');
title(sprintf('Solid acceleration at x = %.2f m', MeshU.coords(Control.plotu,1)));
hold off
if saveGraphs_on
    exportgraphics(gcf,'Vel_time.png','Resolution',300)
end
%% porosity vs depth
if contains(Control.PMmodel, 'UPN')
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
if contains(Control.PMmodel, 'UPN')
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

if contains(Control.PMmodel, 'UPU')
    figure;
    %% fluid displacement vs depth
    subplot(2,2,1);
    plot(MeshU.coords, Solution.uf,'m','LineWidth',2);
    hold on
    xlabel('Column depth [m]');
    ylabel('u (fluid) [m]');
    title(sprintf('Fluid displacement at t = %.0f s', Control.tend));
    hold off
    if saveGraphs_on
        exportgraphics(gcf,'DisplFluid_depth.png','Resolution',300)
    end
    
    %% fluid displacement vs time
    subplot(2,2,3);
    plot(Plot.time, Plot.uf_time,'m','LineWidth',2);
    hold on
    xlabel('Time [s]');
    ylabel('u (fluid) [m]');
    title(sprintf('Fluid displacement at x = %.2f m', MeshU.coords(Control.plotu,1)));
    hold off
    if saveGraphs_on
        exportgraphics(gcf,'DisplFluid_time.png','Resolution',300)
    end
    %% fluid velocity vs time
    subplot(2,2,2);
    plot(Plot.time, Plot.ufdot_time,'c','LineWidth',2);
    hold on
    xlabel('Time [s]');
    ylabel('udot (fluid) [m/s]');
    title(sprintf('Fluid velocity at x = %.2f m', MeshU.coords(Control.plotu,1)));
    hold off
    if saveGraphs_on
        exportgraphics(gcf,'VelFluid_time.png','Resolution',300)
    end
    %% fluid acceleration vs time
    subplot(2,2,4);
    plot(Plot.time, Plot.uf2dot_time,'c','LineWidth',2);
    hold on
    xlabel('Time [s]');
    ylabel('u2dot (fluid) [m/s]');
    title(sprintf('Fluid acceleration at x = %.2f m', MeshU.coords(Control.plotu,1)));
    hold off
    if saveGraphs_on
        exportgraphics(gcf,'VelFluid_time.png','Resolution',300)
    end
end

end