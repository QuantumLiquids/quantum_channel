% Set the values of L, Delta, h, and D
L = 100;
Delta = 1;
h = 0.5;
gamma_values = 0:0.1:1.1;
dump_figure = true;  % Set to true to save the figure as EPS

% Set the directory to search for files
directory = '../data/';
figure_directory = '../note_figure/';

% Initialize arrays to store gamma and corresponding middle bond entanglement entropy values
gamma_values_used = [];
middle_ee_values = [];

for gamma = gamma_values
    % Search for files with the fixed values of L, Delta, h, gamma, and D
    search_pattern = ['ee_brownianL', num2str(L), 'Delta', num2str(Delta), ...
                      'h', num2str(h), 'gamma', num2str(gamma), 'D*'];
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
        tr_parts = strsplit(current_file, 'gamma');
        current_D = str2double(extractAfter(tr_parts{2}, 'D'));
        if current_D > max_D
            max_D = current_D;
            max_D_file = current_file;
        end
    end

    % Load data from the file with the maximum D value for the current gamma
    filename = fullfile(directory, max_D_file);
    data = fread(fopen(filename, 'rb'), 'double');

    % Extract the middle value (middle bond entanglement entropy)
    middle_ee = data(L/2);

    % Store the gamma and middle bond entanglement entropy values in arrays
    gamma_values_used = [gamma_values_used, gamma];
    middle_ee_values = [middle_ee_values, middle_ee];

    % Plot the entanglement entropy
    plot(x, data, 'o'); hold on;
    xlabel('Site','Interpreter','latex')
    ylabel('Entanglement Entropy','Interpreter','latex')
    set(gca,'fontsize',24);
    set(gca,'linewidth',1.5);
    set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
    set(get(gca,'XLabel'),'FontSize',24); 
    set(get(gca,'YLabel'),'FontSize',24);
    % title(['Gamma = ', num2str(gamma)]);
end

if dump_figure
    % Save the figure as EPS file with the specified name
    figure_name_eps = [tr_parts{1}, '.eps'];
    figure_path = fullfile(figure_directory, figure_name_eps);
    saveas(gcf, figure_path, 'epsc');
    disp(['Figure saved as: ', figure_path]);
end

% Plot gamma dependence of middle bond entanglement entropy
figure;
plot(gamma_values_used, middle_ee_values, 'o-', 'LineWidth', 2)
xlabel('$\gamma$','Interpreter','latex')
ylabel('Middle Bond Entanglement Entropy','Interpreter','latex')
set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24);

if dump_figure
    % Save the figure as EPS file with the specified name
    figure_name_eps = 'gamma_dependence_middle_ee.eps';
    figure_path = fullfile(figure_directory, figure_name_eps);
    saveas(gcf, figure_path, 'epsc');
    disp(['Gamma dependence of middle bond entanglement entropy figure saved as: ', figure_path]);
end