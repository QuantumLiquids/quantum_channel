% Set the values of L, Delta, h, and D
L = 100;
Delta = 1;
h = 0.5;
gamma_values = 0:0.1:1.1;
dump_figure = true;  % Set to true to save the figure as EPS

% Set the directory to search for files
directory = '../data/';
figure_directory = '../note_figure/';

% Initialize arrays to store gamma and corresponding decay lengths (1/decay_rate)
gamma_values_used = [];
decay_lengths = [];

for gamma = gamma_values
    % Search for files with the fixed values of L, Delta, h, gamma, and D in filename1
    search_pattern = ['PMbrownianL', num2str(L), 'Delta', num2str(Delta), ...
                      'h', num2str(h), 'gamma', num2str(gamma), 'D*.json'];
    file_list = dir(fullfile(directory, search_pattern));

    if isempty(file_list)
        % Skip this gamma if no suitable data file is found
        disp(['No files found with the specified parameters for gamma = ', num2str(gamma)]);
        continue;
    end

    % Find the file with the maximum D value for the current gamma
    max_D = -Inf;
    max_D_file = '';

    for i = 1:length(file_list)
        current_file = file_list(i).name;
        [~, name, ~] = fileparts(current_file);
        str_parts = strsplit(name, 'D');
        current_D = str2double(str_parts{end});
        if current_D > max_D
            max_D = current_D;
            max_D_file = current_file;
        end
    end

    fprintf("gamma = %.1f, D = %d\n", gamma, max_D);
    % Load data from the file with the maximum D value for the current gamma
    filename1 = fullfile(directory, max_D_file);
    filename2 = strrep(filename1, 'PM', 'MP');  % Generate filename2

    data1 = jsondecode(fileread(filename1));
    data2 = jsondecode(fileread(filename2));
    distances = cell2mat(cellfun(@(x) x{1}(2) - x{1}(1), data1, 'UniformOutput', false));
    u1_values = cell2mat(cellfun(@(x) x{2}, data1, 'UniformOutput', false)) + cell2mat(cellfun(@(x) x{2}, data2, 'UniformOutput', false));

    % Filter the data for distances greater than 20
    filtered_indices = distances > 20;
    filtered_distances = distances(filtered_indices);
    filtered_u1_values = u1_values(filtered_indices);

    % Fit an exponential decay function to the filtered data
    % fit_func = @(p, x) p(1) * exp(-p(2) * x);
    % initial_guess = [1, 0.1]; % Initial guess for the parameters [amplitude, decay rate]
    % fit_result = lsqcurvefit(fit_func, initial_guess, filtered_distances, filtered_u1_values);
    % 
    % decay_rate = fit_result(2);
    % decay_length = 1 / decay_rate;

    % Fit an exponential decay function to the filtered data
    p = fit((filtered_distances),log(filtered_u1_values),'poly1');
    % plot(fit_x, p.p1*fit_x + p.p2,'-.');hold on;

    decay_rate = -p.p1;
    decay_length = 1 / decay_rate;


    % Store the gamma and decay length values in arrays
    gamma_values_used = [gamma_values_used, gamma];
    decay_lengths = [decay_lengths, decay_length];

    % Plot distance dependent of U(1) Order Parameter
    % semilogy(distances, u1_values, 'o')
    % hold on
    % semilogy(filtered_distances, fit_func(fit_result, filtered_distances), 'r', 'LineWidth', 2)
    % hold off
    % xlabel('Distance','Interpreter','latex')
    % ylabel('$U(1)$ Order Parameter','Interpreter','latex')
    % set(gca,'fontsize',24);
    % set(gca,'linewidth',1.5);
    % set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
    % set(get(gca,'XLabel'),'FontSize',24); 
    % set(get(gca,'YLabel'),'FontSize',24);
    % title(['Gamma = ', num2str(gamma)]);

    % if dump_figure
    %     % Save the figure as EPS file with the specified name
    %     figure_name_eps = ['U1brownianL', num2str(L), 'Delta', num2str(Delta), ...
    %                       'h', num2str(h), 'gamma', num2str(gamma), 'D', num2str(max_D), '.eps'];
    %     figure_path = fullfile(figure_directory, figure_name_eps);
    %     saveas(gcf, figure_path, 'epsc');
    %     disp(['Figure saved as: ', figure_path]);
    % end
end

% Plot gamma dependence of decay lengths
figure;
semilogy(gamma_values_used, decay_lengths, 'o-', 'LineWidth', 2)
xlabel('$\gamma$','Interpreter','latex')
ylabel('$U(1)$ order parameter correlation length ($\xi$)','Interpreter','latex')
set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24);

if dump_figure
    % Save the figure as EPS file with the specified name
    figure_name_eps = 'gamma_dependence_u1_correlation_length.eps';
    figure_path = fullfile(figure_directory, figure_name_eps);
    saveas(gcf, figure_path, 'epsc');
    disp(['Gamma dependence of decay length figure saved as: ', figure_path]);
end
