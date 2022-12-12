function PlotFrequencies (phi_u, omega2_u, phi_p, omega2_p, MeshU, MeshP, Control, BC)
% Plot first 6 natural frequencies for u and p

% time
t = Control.t;
% displacement matrix
u = zeros(MeshU.nDOF,6);
% pressure matrix
p = zeros(MeshP.nDOF,6);

% plot first natural frequencies for displacement field
figure;
hold on
sgtitle(sprintf('Displacement Natural Frequencies at t = %.2f s',t ));
for modeu = 1:6
    u(BC.free_u, modeu) = phi_u(:, modeu) * exp(1i*sqrt(omega2_u(modeu,modeu)) * t);
    subplot(2,3, modeu);
    plot(MeshU.coords, u(:,modeu), 'g', 'LineWidth', 2);
    xlabel('Length (m)');
    ylabel('Displacement (m)');
    title(sprintf('Mode %.0f, omega = %.2f Hz', modeu, sqrt(omega2_u(modeu,modeu)) ));
end
hold off

% plot first natural frequencies for pressure field
figure;
hold on
sgtitle(sprintf('Pressure Natural Frequencies at t = %.2f s',t ));
for modep = 1:6
    p(BC.free_p, modep) = phi_p(:, modep) * exp(1i*sqrt(omega2_p(modep,modep)) * t);
    subplot(2,3, modep);
    plot(MeshP.coords, p(:,modep), 'b', 'LineWidth', 2);
    xlabel('Length (m)');
    ylabel('Pressure (Pa)');
    title(sprintf('Mode %.0f, omega = %.2f Hz', modep, sqrt(omega2_p(modep,modep)) ));
end
hold off

% plot first mode shapes for displacement field
figure;
hold on
sgtitle(sprintf('Displacement mode shapes at t = %.2f s',t ));
for modep = 1:6
    p(BC.free_p, modep) = phi_p(:, modep) * exp(1i*sqrt(omega2_p(modep,modep)) * t);
    subplot(2,3, modep);
    plot(MeshP.coords, p(:,modep), 'b', 'LineWidth', 2);
    xlabel('Length (m)');
    ylabel('Pressure (Pa)');
    title(sprintf('Mode %.0f, omega = %.2f Hz', modep, sqrt(omega2_p(modep,modep)) ));
end
hold off


end