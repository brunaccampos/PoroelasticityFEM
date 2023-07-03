function PlotSynthetics(MeshU, MeshP, MeshN, Plot, Control)
% ------------------------------------------------------------------------
% Plot synthetics
% ------------------------------------------------------------------------

% solid displacement
figure;
for i = 1:length(Control.ploturow)
    plot(Plot.time, Plot.u_synthetic(:,i) + MeshU.coords(Control.ploturow(i)./2,2),'k','LineWidth',1.5);
    hold on
end
xlabel('Time (s)');
ylabel('Coordinate in y (m)');
title(sprintf('Synthetics Solid displacement at x = %.2f m', Control.depthplot));
hold off

% fluid pressure
figure;
for i = 1:length(Control.plotprow)
    plot(Plot.time, Plot.p_synthetic(:,i) + MeshP.coords(Control.plotprow(i),2),'k','LineWidth',1.5);
    hold on
end
xlabel('Time (s)');
ylabel('Coordinate in y (m)');
title(sprintf('Synthetics Fluid pressure at x = %.2f m', Control.depthplot));
hold off

% porosity
if contains(Control.PMmodel, 'UPN')
    figure;
    for i = 1:length(Control.plotprow)
        plot(Plot.time, Plot.n_synthetic(:,i) + MeshN.coords(Control.plotprow(i),2),'k','LineWidth',1.5);
        hold on
    end
    xlabel('Time (s)');
    ylabel('Coordinate in y (m)');
    title(sprintf('Synthetics Porosity at x = %.2f m', Control.depthplot));
    hold off
end

% fluid displacement
if contains(Control.PMmodel, 'UPU')
    figure;
    for i = 1:length(Control.ploturow)
        plot(Plot.time, Plot.uf_synthetic(:,i) + MeshU.coords(Control.ploturow(i)./2,2),'k','LineWidth',1.5);
        hold on
    end
    xlabel('Time (s)');
    ylabel('Coordinate in y (m)');
    title(sprintf('Synthetics Fluid displacement at x = %.2f m', Control.depthplot));
    hold off
end

end