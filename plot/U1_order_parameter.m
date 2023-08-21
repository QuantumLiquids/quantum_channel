% Set the values of L, Delta, h, gamma, and D
L = 100;
Delta = 1;
h = 0.5;
gamma = 0.1;
D = 1500;
dump_figure = true;  % Set to true to save the figure as EPS

filename1 = ['../data/PMbrownianL', num2str(L), 'Delta', num2str(Delta), ...
            'h', num2str(h), 'gamma', num2str(gamma), 'D', num2str(D),'.json'];

filename2 = ['../data/PMbrownianL', num2str(L), 'Delta', num2str(Delta), ...
            'h', num2str(h), 'gamma', num2str(gamma), 'D', num2str(D),'.json'];

% Load data from file
data1 = jsondecode(fileread(filename1));
data2 = jsondecode(fileread(filename2));
distances = cell2mat(cellfun(@(x) x{1}(2) - x{1}(1), data1, 'UniformOutput', false));
u1_values = cell2mat(cellfun(@(x) x{2}, data1, 'UniformOutput', false)) + cell2mat(cellfun(@(x) x{2}, data2, 'UniformOutput', false));

% Filter the data for distances greater than 20
filtered_indices = distances > 20;
filtered_distances = distances(filtered_indices);
filtered_u1_values = u1_values(filtered_indices);

% Fit an exponential decay function to the filtered data
fit_func = @(p, x) (p(1) * exp(-p(2) * x)  );
initial_guess = [1,0.1,0]; % Initial guess for the parameters [amplitude, decay rate]
fit_result = lsqcurvefit(fit_func, initial_guess, filtered_distances, filtered_u1_values);

amplitude = fit_result(1);
decay_rate = fit_result(2);
% infty_correlation = fit_result(3);

% Plot distance dependent of S2_order_parameter
semilogy(distances, u1_values, 'o')
hold on
semilogy(filtered_distances, fit_func(fit_result, filtered_distances), 'r', 'LineWidth', 2)
% hold off
xlabel('Distance','Interpreter','latex')
ylabel('$U(1)$ Order Parameter','Interpreter','latex')
set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24);

% Print parameter values and decay length in the figure
text(0.7, 0.5, ['$\Delta =$ ', num2str(Delta)], 'Interpreter', 'latex', 'FontSize', 24, 'Units', 'normalized')
text(0.7, 0.6, ['$h =$ ', num2str(h)], 'Interpreter', 'latex', 'FontSize', 24, 'Units', 'normalized')
text(0.7, 0.7, ['$\gamma =$ ', num2str(gamma)], 'Interpreter', 'latex', 'FontSize', 24, 'Units', 'normalized')
text(0.7, 0.8, ['$\xi =$ ', num2str(1/decay_rate)], 'Interpreter', 'latex', 'FontSize', 24, 'Units', 'normalized')
% text(0.7, 0.9, ['correlation($\infty)', num2str(infty_correlation) , 'Interpreter', 'latex', 'FontSize', 24, 'Units', 'normalized'])


if dump_figure
    % Save the figure as EPS file with the specified name
    figure_name_eps = ['U1brownianL', num2str(L), 'Delta', num2str(Delta), ...
                      'h', num2str(h), 'gamma', num2str(gamma), 'D', num2str(D), '.eps'];
    figure_path = fullfile('../note_figure', figure_name_eps);
    saveas(gcf, figure_path, 'epsc');
    disp(['Figure saved as: ', figure_path]);
end
