% Set the values of L, Delta, h, and D
L = 100;
Delta = 1;
h = 0.5;
gamma_values = sort([0:0.1:1.1, 0.82:0.03:0.97]);
dump_figure = true;  % Set to true to save the figure as EPS

% Set the directory to search for files
directory = '../data/';
figure_directory = '../note_figure/';

% Initialize arrays to store gamma and corresponding xi values
gamma_values_used = [];
order_values = [];

for gamma = gamma_values
    % Search for files with the fixed values of L, Delta, h, gamma, and D
    search_pattern = ['S2symbrownianL', num2str(L), 'Delta', num2str(Delta), ...
        'h', num2str(h), 'gamma', num2str(gamma), 'D*.json'];
    file_list = dir(fullfile(directory, search_pattern));

    if isempty(file_list)
        % Skip this gamma if no suitable data file is found
        disp(['No files found with the specified parameters for gamma = ', num2str(gamma)]);
        continue;
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
    s2_values = cell2mat(cellfun(@(x) x{2}, data, 'UniformOutput', false));
    order = mean(s2_values);
    % Store the gamma and xi values in arrays
    gamma_values_used = [gamma_values_used, gamma];
    order_values = [order_values, order];


    fprintf('gamma = %.2f:, Max D = %d, order_value = %.5f\n ', (gamma), max_D, order);

end

% Plot gamma dependence of xi
% figure;
plot(gamma_values_used, order_values, 'o-', 'LineWidth', 2)
xlabel('$\gamma$','Interpreter','latex')
ylabel('$S^2$','Interpreter','latex')
set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24);

if dump_figure
    % Save the figure as EPS file with the specified name
    figure_name_eps = 'gamma_dependence_xi.eps';
    figure_path = fullfile(figure_directory, figure_name_eps);
    saveas(gcf, figure_path, 'epsc');
    disp(['Gamma dependence of xi figure saved as: ', figure_path]);
end
