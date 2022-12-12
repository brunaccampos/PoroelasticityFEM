function postprocessor(Mesh, Control, Plot, u, p)
% Post processor

if Control.steady
    % displacement vs depth plot
    figure;
    plot(Mesh.coords_u, u,'b','LineWidth',2);
    hold on
    plot(Mesh.coords_u, Plot.u_an,'m--','LineWidth',2);
    xlabel('Column depth (m)');
    ylabel('u (m)');
    title(sprintf('Displacement at t = %.0f s', Control.tend));
    legend('Numerical', 'Analytical');
    hold off
%     exportgraphics(gcf,'Displacement.png','Resolution',300)

    % pressure vs depth plot
    figure;
    plot(Mesh.coords_p, p,'b','LineWidth',2);
    hold on
    plot(Mesh.coords_p, Plot.p_an,'m--','LineWidth',2);
    xlabel('Column depth (m)');
    ylabel('p (Pa)');
    title(sprintf('Pressure at t = %.0f s', Control.tend));
    legend('Numerical', 'Analytical');
    hold off
%     exportgraphics(gcf,'Pressure.png','Resolution',300)

    % pressure vs time plot
    figure;
    plot(Plot.time, Plot.p,'b','LineWidth',2);
    hold on
    plot(Plot.time, Plot.pan,'m--','LineWidth',2);
    xlabel('Time (s)');
    ylabel('p (Pa)');
    title(sprintf('Pressure at x = %.2f m', Mesh.coords_p(1, Control.plotp)));
    legend('Numerical', 'Analytical');
    hold off

else
    % displacement vs depth plot
    figure;
%     subplot(2,1,1);
    plot(Mesh.coords_u,u,'b','LineWidth',2);
    hold on
    xlabel('Column depth (m)');
    ylabel('u (m)');
    title(sprintf('Displacement at t = %.0f s', Control.tend));
    hold off
%     exportgraphics(gcf,'Displacement.png','Resolution',300)

    % pressure vs depth plot
%     subplot(2,1,2);
    figure;
    plot(Mesh.coords_p,p,'b','LineWidth',2);
    hold on
    xlabel('Column depth (m)');
    ylabel('p (Pa)');
    title(sprintf('Pressure at t = %.0f s', Control.tend));
    hold off
%     exportgraphics(gcf,'Pressure.png','Resolution',300)

    % pressure vs time plot
    figure;
%     subplot(2,2,1);
    plot(Plot.time, Plot.p,'b','LineWidth',2);
    hold on
    xlabel('Time (s)');
    ylabel('p (Pa)');
    title(sprintf('Pressure at x = %.2f m', Mesh.coords_p(1, Control.plotp)));
    hold off
%     exportgraphics(gcf,'PressureTime.png','Resolution',300)

    % velocity vs time plot
%     subplot(2,2,2);
    figure;
    plot(Plot.time, Plot.udot,'b','LineWidth',2);
    hold on
    xlabel('Time (s)');
    ylabel('udot (m/s)');
    title(sprintf('Solid skeleton velocity at x = %.2f m', Mesh.coords_u(1, Control.plotu)));
    hold off
%     exportgraphics(gcf,'VelocityTime.png','Resolution',300)

    % displacement vs time plot
%     subplot(2,2,3);
    figure;
    plot(Plot.time, Plot.u,'b','LineWidth',2);
    hold on
    xlabel('Time (s)');
    ylabel('u (m/s)');
    title(sprintf('Solid skeleton displacement at x = %.2f m', Mesh.coords_u(1, Control.plotu)));
    hold off
%     exportgraphics(gcf,'DisplacementTime.png','Resolution',300)
end

end