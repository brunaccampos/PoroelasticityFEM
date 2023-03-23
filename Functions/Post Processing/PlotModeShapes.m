function PlotModeShapes(phi_u, omega2_u, phi_p, omega2_p, MeshU, MeshP, Control, BC, config_name, vtk_dir)
% Plot first 6 natural frequencies for u and p
% ------------------------------------------------------------------------

% sort modes first
omega2_u_vec = diag(omega2_u);
[omega2_u_sorted, index] = sort(omega2_u_vec);
phi_u_sorted = zeros(length(phi_u), length(phi_u));
for i = 1: length(phi_u)
    phi_u_sorted(:,i) = phi_u(:,index(i));
end

omega2_p_vec = diag(omega2_p);
[omega2_p_sorted, index] = sort(omega2_p_vec);
phi_p_sorted = zeros(length(phi_p), length(phi_p));
for i = 1: length(phi_p)
    phi_p_sorted(:,i) = phi_p(:,index(i));
end

if MeshU.nsd == 1
    if Control.steady
        % 1D transient case
        PlotModes1DTransient(phi_u_sorted, omega2_u_sorted, phi_p_sorted, omega2_p_sorted, MeshU, MeshP, Control, BC);
    else
        % 1D dynamic case
        PlotModes1DDynamic(phi_u_sorted, omega2_u_sorted, phi_p_sorted, omega2_p_sorted, MeshU, MeshP, Control, BC);
    end
else
    plotmodes = 6;
    for mode = 1:plotmodes
        PostProcessingModes(mode, phi_u_sorted, phi_p_sorted, MeshU, MeshP, BC, config_name, vtk_dir);
    end
end
end

%% 1D transient plot
function PlotModes1DTransient(phi_u, omega2_u, phi_p, omega2_p, MeshU, MeshP, Control, BC)
% plot first natural mode shapes
figure;
hold on
sgtitle(sprintf('PM Mode Shapes at t = %.2f s', Control.tend));
modeshapeu = zeros(MeshU.nDOF,6);
modeshapep = zeros(MeshP.nDOF,6);
for mode = 1:6
    modeshapeu(BC.free_u, mode) = phi_u(:,mode);
    modeshapep(BC.free_p, mode) = phi_p(:,mode);
    subplot(2,3, mode);
    plot(MeshU.coords, modeshapeu(:,mode), 'b', 'LineWidth', 2);
    hold on
    plot(MeshP.coords, modeshapep(:,mode), 'k', 'LineWidth', 2);
    xlabel('Length (m)');
    ylabel('Shape');
    legend('u', 'p');
    title(sprintf('# %.0f, omU = %.2f Hz, omP = %.2d Hz', mode, sqrt(omega2_u(mode)*1e9), sqrt(omega2_p(mode)*1e-9)));
end
hold off
end

%% 1D dynamic plot
function PlotModes1DDynamic(phi_u, omega2_u, phi_p, omega2_p, MeshU, MeshP, Control, BC)
% plot first natural mode shapes
figure;
hold on
sgtitle(sprintf('PM Mode Shapes at t = %.2f s', Control.tend));
modeshapeu = zeros(MeshU.nDOF,6);
modeshapep = zeros(MeshP.nDOF,6);
for mode = 1:6
    modeshapeu(BC.free_u, mode) = phi_u(:,mode) / norm(phi_u(:,mode));
    modeshapep(BC.free_p, mode) = phi_p(:,mode);
    subplot(2,3, mode);
    plot(MeshU.coords, modeshapeu(:,mode), 'b', 'LineWidth', 2);
    hold on
    plot(MeshP.coords, modeshapep(:,mode), 'k', 'LineWidth', 2);
    xlabel('Length (m)');
    ylabel('Shape');
    legend('u', 'p');
    title(sprintf('# %.0f, omU = %.2f Hz, omP = %.2d Hz', mode, sqrt(omega2_u(mode)), sqrt(omega2_p(mode)*1e-9)));
end
hold off

end
