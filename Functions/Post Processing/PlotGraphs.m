function plotGraphs(MeshU, MeshP, Control, Plot, u, p, saveGraphs_on)
% Plot results in graphs
% Input: Mesh, Control, Plot, u, p
% Output: plots of
%               displacement vs depth
%               pressure vs depth
%               displacement vs time at specific node
%               pressure vs time at specific node

%% Quasi-steady state case
if Control.steady
    % displacement vs depth plot
    figure;
    plot(MeshU.coords, u,'b','LineWidth',2);
    hold on
    plot(MeshU.coords, Plot.u_an,'m--','LineWidth',2);
    xlabel('Column depth (m)');
    ylabel('u (m)');
    title(sprintf('Displacement at t = %.0f s', Control.tend));
    legend('Numerical', 'Analytical');
    hold off
    if saveGraphs_on
        exportgraphics(gcf,'Displ_depth.png','Resolution',300)
    end

    % pressure vs depth plot
    figure;
    plot(MeshP.coords, p,'b','LineWidth',2);
    hold on
    plot(MeshP.coords, Plot.p_an,'m--','LineWidth',2);
    xlabel('Column depth (m)');
    ylabel('p (Pa)');
    title(sprintf('Pressure at t = %.0f s', Control.tend));
    legend('Numerical', 'Analytical');
    hold off
    if saveGraphs_on
        exportgraphics(gcf,'Press_depth.png','Resolution',300)
    end

    % pressure vs time plot
    figure;
    plot(Plot.time, Plot.p,'b','LineWidth',2);
    hold on
    plot(Plot.time, Plot.pan,'m--','LineWidth',2);
    xlabel('Time (s)');
    ylabel('p (Pa)');
    title(sprintf('Pressure at x = %.2f m', MeshP.coords(Control.plotp,1)));
    legend('Numerical', 'Analytical');
    hold off
    if saveGraphs_on
        exportgraphics(gcf,'Press_time.png','Resolution',300)
    end

    %% Transient case
else
    % 1D case
    if MeshU.nsd == 1
        % displacement vs depth plot
        figure;
        plot(MeshU.coords,u,'b','LineWidth',2);
        hold on
        xlabel('Column depth (m)');
        ylabel('u (m)');
        title(sprintf('Displacement at t = %.0f s', Control.tend));
        hold off
        if saveGraphs_on
            exportgraphics(gcf,'Displ_depth.png','Resolution',300)
        end

        % pressure vs depth plot
        figure;
        plot(MeshP.coords,p,'b','LineWidth',2);
        hold on
        xlabel('Column depth (m)');
        ylabel('p (Pa)');
        title(sprintf('Pressure at t = %.0f s', Control.tend));
        hold off
        if saveGraphs_on
            exportgraphics(gcf,'Press_depth.png','Resolution',300)
        end

        % pressure vs time plot
        figure;
        plot(Plot.time, Plot.p,'b','LineWidth',2);
        hold on
        xlabel('Time (s)');
        ylabel('p (Pa)');
        title(sprintf('Pressure at x = %.2f m', MeshP.coords(Control.plotp,1)));
        hold off
        if saveGraphs_on
            exportgraphics(gcf,'Press_time.png','Resolution',300)
        end

        % velocity vs time plot
        figure;
        plot(Plot.time, Plot.udot,'b','LineWidth',2);
        hold on
        xlabel('Time (s)');
        ylabel('udot (m/s)');
        title(sprintf('Solid skeleton velocity at x = %.2f m', MeshU.coords(Control.plotu,1)));
        hold off
        if saveGraphs_on
            exportgraphics(gcf,'Vel_time.png','Resolution',300)
        end

        % displacement vs time plot
        figure;
        plot(Plot.time, Plot.u,'b','LineWidth',2);
        hold on
        xlabel('Time (s)');
        ylabel('u (m/s)');
        title(sprintf('Solid skeleton displacement at x = %.2f m', MeshU.coords(Control.plotu,1)));
        hold off
        if saveGraphs_on
            exportgraphics(gcf,'Displ_time.png','Resolution',300)
        end

        % 2D case
    else
        % displacement vs depth plot
        figure;
        scatter(MeshU.coords(:,1), MeshU.coords(:,2),'c','filled');
        hold on
        scatter(MeshU.coords(:,1)+u(1:2:end), MeshU.coords(:,2)+u(2:2:end),'m', 'filled');
        xlabel('x (m)');
        ylabel('y (m)');
        title(sprintf('Displacement at t = %.0f s', Control.tend));
        hold off
        if saveGraphs_on
            exportgraphics(gcf,'Displ_depth.png','Resolution',300)
        end

        % pressure vs time plot
        figure;
        plot(Plot.time, Plot.p,'b','LineWidth',2);
        hold on
        xlabel('Time (s)');
        ylabel('p (Pa)');
        title(sprintf('Pressure at x = %.2f m, y = %.2f m', MeshP.coords(Control.plotp,1), MeshP.coords(Control.plotp,2)));
        hold off
        if saveGraphs_on
            exportgraphics(gcf,'Press_time.png','Resolution',300)
        end

        % velocity vs time plot
        figure;
        plot(Plot.time, Plot.udot,'b','LineWidth',2);
        hold on
        xlabel('Time (s)');
        ylabel('udot (m/s)');
        title(sprintf('Solid skeleton velocity at x = %.2f m, y = %.2f m', MeshU.coords(Control.plotu/2,1), MeshU.coords(Control.plotu/2,2)));
        hold off
        if saveGraphs_on
            exportgraphics(gcf,'Vel_time.png','Resolution',300)
        end

        % displacement vs time plot
        figure;
        plot(Plot.time, Plot.u,'b','LineWidth',2);
        hold on
        xlabel('Time (s)');
        ylabel('u (m/s)');
        title(sprintf('Solid skeleton displacement at x = %.2f m, y = %.2f m', MeshU.coords(Control.plotu/2,1), MeshU.coords(Control.plotu/2,2)));
        hold off
        if saveGraphs_on
            exportgraphics(gcf,'Displ_time.png','Resolution',300)
        end
    end


end

end