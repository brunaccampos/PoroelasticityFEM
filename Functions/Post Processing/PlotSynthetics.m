function PlotSynthetics(MeshU, MeshP, MeshN, Plot, Control)
% ------------------------------------------------------------------------
% Plot synthetics at fixed coordinate (x or y) in 2D domain
% ------------------------------------------------------------------------

% initialize figure
figure;
tiledlayout(2,3);

%% solid displacement
nexttile
for i = 1:length(Control.ploturow)
    plot(Plot.u_synthetic(:,i) + MeshU.coords(Control.ploturow(i)./Control.DOFplot, Control.depthDir), Plot.time, 'k','LineWidth',1.5);
    hold on
end
ylabel('Time [s]');
xlabel('Coordinate [m]');
xlim([0 max(MeshU.coords(:,Control.depthDir))]);
title(sprintf('Solid displ. at fixed %.2f m, DOF %.0f', Control.depthplot, Control.DOFplot));
hold off

%% solid velocity
% normalized values
udot_normalized = Plot.udot_synthetic./max(Plot.udot_synthetic,[],'all');

nexttile
for i = 1:length(Control.ploturow)
    plot(udot_normalized(:,i) + MeshU.coords(Control.ploturow(i)./Control.DOFplot, Control.depthDir), Plot.time, 'k','LineWidth',1.5);
    hold on
end
ylabel('Time [s]');
xlabel('Coordinate [m]');
xlim([0 max(MeshU.coords(:,Control.depthDir))]);
title(sprintf('Solid vel. at fixed %.2f m, DOF %.0f', Control.depthplot, Control.DOFplot));
hold off

%% fluid pressure
% normalized values
p_normalized = Plot.p_synthetic./max(abs(Plot.p_synthetic),[],'all');

nexttile
for i = 1:length(Control.plotprow)
    plot(p_normalized(:,i) + MeshP.coords(Control.plotprow(i), Control.depthDir), Plot.time, 'k','LineWidth',1.5);
    hold on
end
ylabel('Time [s]');
xlabel('Coordinate [m]');
xlim([0 max(MeshU.coords(:,Control.depthDir))]);
title(sprintf('Fluid press. at fixed %.2f m', Control.depthplot));
hold off

%% porosity
if contains(Control.PMmodel, 'UPN')
    nexttile
    for i = 1:length(Control.plotprow)
        plot(Plot.n_synthetic(:,i) + MeshN.coords(Control.plotprow(i), Control.depthDir), Plot.time, 'k','LineWidth',1.5);
        hold on
    end
    ylabel('Time [s]');
    xlabel('Coordinate [m]');
    xlim([0 max(MeshU.coords(:,Control.depthDir))]);
    title(sprintf('Porosity at fixed %.2f m, DOF %.0f', Control.depthplot, Control.DOFplot));
    hold off
end

if contains(Control.PMmodel, 'UPU')
    %% fluid displacement
    nexttile
    for i = 1:length(Control.ploturow)
        plot(Plot.uf_synthetic(:,i) + MeshU.coords(Control.ploturow(i)./Control.DOFplot, Control.depthDir), Plot.time, 'k','LineWidth',1.5);
        hold on
    end
    ylabel('Time [s]');
    xlabel('Coordinate [m]');
    xlim([0 max(MeshU.coords(:,Control.depthDir))]);
    title(sprintf('Fluid displ. at fixed %.2f m, DOF %.0f', Control.depthplot, Control.DOFplot));
    hold off

    %% fluid velocity
    % normalized values
    ufdot_normalized = Plot.ufdot_synthetic./max(Plot.ufdot_synthetic,[],'all');

    nexttile
    for i = 1:length(Control.ploturow)
        plot(ufdot_normalized(:,i) + MeshU.coords(Control.ploturow(i)./Control.DOFplot, Control.depthDir), Plot.time, 'k','LineWidth',1.5);
        hold on
    end
    ylabel('Time [s]');
    xlabel('Coordinate [m]');
    xlim([0 max(MeshU.coords(:,Control.depthDir))]);
    title(sprintf('Fluid vel. at fixed %.2f m, DOF %.0f', Control.depthplot, Control.DOFplot));
    hold off
end

end