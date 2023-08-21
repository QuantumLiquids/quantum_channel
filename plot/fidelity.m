% Parameters
L_values = [32,64, 100];
Delta = 1;
h = 0.5;

% Initialize the figure
figure;

% Loop over different L values
for i = 1:length(L_values)
    L = L_values(i);
    filename = ['../data/fidelityL', num2str(L), 'Delta', num2str(Delta), ...
                          'h', num2str(h), '.csv'];
    data = csvread(filename, 1, 0); % Skip the header (first row) while reading

    % Extract L, Delta, h from the first row
    L = data(1, 1);
    Delta = data(1, 2);
    h = data(1, 3);

    % Extract the gamma pairs and fidelity values (starting from the second row)
    gamma1 = data(2:end, 1);
    gamma2 = data(2:end, 2);
    fidelity_data = data(2:end, 3);

    % Plot the curve
    plot((gamma1 + gamma2) / 2, fidelity_data, '-o', 'LineWidth', 2);
    hold on;
end

% Add legend and labels
l = legend(cellstr(num2str(L_values', '$L = %d$')), 'Interpreter', 'latex');
set(l,'Box','off');
set(l,'Fontsize',24);
set(l,'Location','SouthWest');

xlabel('$\gamma$', 'Interpreter', 'latex')
ylabel('fidelity', 'Interpreter', 'latex')
set(gca, 'fontsize', 24);
set(gca, 'linewidth', 1.5);
set(get(gca, 'Children'), 'linewidth', 2); % Set line width to 2 points
set(get(gca, 'XLabel'), 'FontSize', 24); 
set(get(gca, 'YLabel'), 'FontSize', 24);

% Save the figure as fidelity.eps
saveas(gcf, '../note_figure/fidelity.eps', 'epsc');
