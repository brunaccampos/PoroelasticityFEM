function PlotSynthetics(MeshU, MeshP, MeshN, Plot, Control)
% ------------------------------------------------------------------------
% Plot synthetics at fixed coordinate (x or y) in 2D domain
% ------------------------------------------------------------------------

% scaling factor
scale = 1;

% initialize figure
figure;
tiledlayout(3,3);

%% solid displacement in x
% find half of array
half = ceil(length(Control.ploturow)/2);
% normalized values
ux_normalized = scale * Plot.u_synthetic(:,1:half)./max(abs(Plot.u_synthetic(:,1:half)),[],'all');

nexttile
for i = 1:half
    plot(ux_normalized(:,i) + MeshU.coords((Control.ploturow(i)+1)./2, Control.depthDir), Plot.time, 'k','LineWidth',1.5);
    hold on
end
ylabel('Time [s]');
xlabel('Coordinate [m]');
% xlim([0 max(MeshU.coords(:,Control.depthDir))]);
title(sprintf('Solid displ. in x at fixed %.0f, %.2f m', Control.depthDir, Control.depthplot));
hold off

%% solid velocity in x
% normalized values
uxdot_normalized = scale * Plot.udot_synthetic(:,1:half)./max(abs(Plot.udot_synthetic(:,1:half)),[],'all');

nexttile
for i = 1:half
    plot(uxdot_normalized(:,i) + MeshU.coords((Control.ploturow(i)+1)./2, Control.depthDir), Plot.time, 'k','LineWidth',1.5);
    hold on
end
ylabel('Time [s]');
xlabel('Coordinate [m]');
% xlim([0 max(MeshU.coords(:,Control.depthDir))]);
title(sprintf('Solid vel. in x at fixed %.0f, %.2f m', Control.depthDir, Control.depthplot));
hold off

%% fluid pressure
% normalized values
p_normalized = scale * Plot.p_synthetic./max(abs(Plot.p_synthetic),[],'all');

nexttile
for i = 1:length(Control.plotprow)
    plot(p_normalized(:,i) + MeshP.coords(Control.plotprow(i), Control.depthDir), Plot.time, 'k','LineWidth',1.5);
    hold on
end
ylabel('Time [s]');
xlabel('Coordinate [m]');
% xlim([0 max(MeshU.coords(:,Control.depthDir))]);
title(sprintf('Fluid press. at fixed %.0f, %.2f m', Control.depthDir, Control.depthplot));
hold off

%% solid displacement in y
% normalized values
uy_normalized = scale * Plot.u_synthetic(:,half+1:end)./max(abs(Plot.u_synthetic(:,half+1:end)),[],'all');

nexttile
for i = 1:half
    plot(uy_normalized(:,i) + MeshU.coords((Control.ploturow(i+half))./2, Control.depthDir), Plot.time, 'k','LineWidth',1.5);
    hold on
end
ylabel('Time [s]');
xlabel('Coordinate [m]');
% xlim([0 max(MeshU.coords(:,Control.depthDir))]);
title(sprintf('Solid displ. in y at fixed %.0f, %.2f m', Control.depthDir, Control.depthplot));
hold off

%% solid velocity in y
% normalized values
uydot_normalized = scale * Plot.udot_synthetic(:,half+1:end)./max(abs(Plot.udot_synthetic(:,half+1:end)),[],'all');

nexttile
for i = 1:half
    plot(uydot_normalized(:,i) + MeshU.coords((Control.ploturow(i+half))./2, Control.depthDir), Plot.time, 'k','LineWidth',1.5);
    hold on
end
ylabel('Time [s]');
xlabel('Coordinate [m]');
% xlim([0 max(MeshU.coords(:,Control.depthDir))]);
title(sprintf('Solid vel. in y at fixed %.0f, %.2f m', Control.depthDir, Control.depthplot));
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
    title(sprintf('Porosity at fixed %.0f, %.2f m', Control.depthDir, Control.depthplot));
    hold off
end

if contains(Control.PMmodel, 'UPU')
    %% fluid displacement in x
    % normalized values
    ufx_normalized = scale * Plot.uf_synthetic(:,1:half)./max(abs(Plot.uf_synthetic(:,1:half)),[],'all');
    
    nexttile
    for i = 1:half
        plot(ufx_normalized(:,i) + MeshU.coords((Control.ploturow(i)+1)./2, Control.depthDir), Plot.time, 'k','LineWidth',1.5);
        hold on
    end
    ylabel('Time [s]');
    xlabel('Coordinate [m]');
    xlim([0 max(MeshU.coords(:,Control.depthDir))]);
    title(sprintf('Fluid displ. in x at fixed %.0f, %.2f m', Control.depthDir, Control.depthplot));
    hold off

    %% fluid velocity in x
    % normalized values
    ufxdot_normalized = scale * Plot.ufdot_synthetic(:,1:half)./max(abs(Plot.ufdot_synthetic(:,1:half)),[],'all');

    nexttile
    for i = 1:half
        plot(ufxdot_normalized(:,i) + MeshU.coords((Control.ploturow(i)+1)./2, Control.depthDir), Plot.time, 'k','LineWidth',1.5);
        hold on
    end
    ylabel('Time [s]');
    xlabel('Coordinate [m]');
    xlim([0 max(MeshU.coords(:,Control.depthDir))]);
    title(sprintf('Fluid vel. in x at fixed %.0f, %.2f m', Control.depthDir, Control.depthplot));
    hold off
    
    %% fluid displacement in y
    % normalized values
    ufy_normalized = scale * Plot.uf_synthetic(:,half+1:end)./max(abs(Plot.uf_synthetic(:,half+1:end)),[],'all');
    
    nexttile
    for i = 1:half
        plot(ufy_normalized(:,i) + MeshU.coords((Control.ploturow(i+half))./2, Control.depthDir), Plot.time, 'k','LineWidth',1.5);
        hold on
    end
    ylabel('Time [s]');
    xlabel('Coordinate [m]');
    xlim([0 max(MeshU.coords(:,Control.depthDir))]);
    title(sprintf('Fluid displ. in y at fixed %.0f, %.2f m', Control.depthDir, Control.depthplot));
    hold off

    %% fluid velocity in y
    % normalized values
    ufydot_normalized = scale * Plot.ufdot_synthetic(:,half+1:end)./max(abs(Plot.ufdot_synthetic(:,half+1:end)),[],'all');
    
    nexttile
    for i = 1:half
        plot(ufydot_normalized(:,i) + MeshU.coords((Control.ploturow(i+half))./2, Control.depthDir), Plot.time, 'k','LineWidth',1.5);
        hold on
    end
    ylabel('Time [s]');
    xlabel('Coordinate [m]');
    xlim([0 max(MeshU.coords(:,Control.depthDir))]);
    title(sprintf('Fluid vel. in y at fixed %.0f, %.2f m', Control.depthDir, Control.depthplot));
    hold off
end

end