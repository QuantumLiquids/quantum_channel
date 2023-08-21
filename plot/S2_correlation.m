% Set the values of L, Delta, h, gamma, and D
L = 100;
Delta = 1;
h = 0.5;
gamma = 0.1;
dump_figure = true;  % Set to true to save the figure as EPS

% Set the directory to search for files
directory = '../data/';
figure_directory = '../note_figure/';

% Search for files with the fixed values of L, Delta, h, gamma
search_pattern = ['S2symbrownianL', num2str(L), 'Delta', num2str(Delta), ...
    'h', num2str(h), 'gamma', num2str(gamma), 'D*.json'];
file_list = dir(fullfile(directory, search_pattern));

if isempty(file_list)
    error('No files found with the specified parameters.')
end

% Find the file with the maximum D value
max_D = -Inf;
max_D_file = '';

for i = 1:length(file_list)
    current_file = file_list(i).name;
    [~, name, ~] = fileparts(current_file);
    str_parts = strsplit(name, 'gamma');
    current_D = str2double(extractAfter(str_parts{2}, 'D'));
    if current_D > max_D
        max_D = current_D;
        max_D_file = current_file;
    end
end

% Load data from the file with the maximum D value
filename = fullfile(directory, max_D_file);
data = jsondecode(fileread(filename));
distances = cell2mat(cellfun(@(x) x{1}(2) - x{1}(1), data, 'UniformOutput', false));
s2_values = cell2mat(cellfun(@(x) x{2}, data, 'UniformOutput', false));

% Filter the data for distances greater than 20
filtered_indices = distances > 20;
filtered_distances = distances(filtered_indices);
filtered_s2_values = s2_values(filtered_indices);

% Fit an exponential decay function to the filtered data
fit_func = @(p, x) p(1) * exp(-p(2) * x);
initial_guess = [1, 0.1]; % Initial guess for the parameters [amplitude, decay rate]
fit_result = lsqcurvefit(fit_func, initial_guess, filtered_distances, filtered_s2_values);

amplitude = fit_result(1);
decay_rate = fit_result(2);

% Plot distance dependent of S2_order_parameter
semilogy(distances, s2_values, 'o')
hold on
semilogy(filtered_distances, fit_func(fit_result, filtered_distances), 'r', 'LineWidth', 2)
hold off

xlabel('Distance','Interpreter','latex')
ylabel('$S^2$ Order Parameter','Interpreter','latex')
set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24);

text(0.7, 0.8, ['$\xi =$ ', num2str(1/decay_rate)], 'Interpreter', 'latex', 'FontSize', 24, 'Units', 'normalized')


% Output the maximum D value
disp(['Maximum D value: ', num2str(max_D)]);

if dump_figure
    % Save the figure as EPS file with max_D_file as the name
    [~, figure_name, ~] = fileparts(max_D_file);
    figure_name_eps = [figure_name, '.eps'];
    figure_path = fullfile(figure_directory, figure_name_eps);
    saveas(gcf, figure_path, 'epsc');
    disp(['Figure saved as: ', figure_path]);
end
