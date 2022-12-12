function Plot2csv (config_name, vtk_dir, Plot)

% file name
filename = fullfile(vtk_dir, [config_name '.xls']);

% step
step = (1:length(Plot.time)).';
% time
time = Plot.time.';
% solid displacement
u = Plot.u.';
% solid velocity
udot = Plot.udot.';
% fluid pressure
p = Plot.p.';

Results = [step, time, u, udot, p];

writematrix(Results, filename);

end