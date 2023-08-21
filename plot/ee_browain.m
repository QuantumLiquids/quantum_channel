% Set the values of L, Delta, h, gamma, and D
L = 100;
Delta = 1;
h = 0.5;
gamma = 0.7;
D = 2500;
dump_figure = true;  % Set to true to save the figure as EPS

% Construct the filename based on the variable values
filename = ['../data/ee_brownianL', num2str(L), 'Delta', num2str(Delta), ...
    'h', num2str(h), 'gamma', num2str(gamma), 'D', num2str(D)];

% Load data from the file
data = fread(fopen(filename, 'rb'), 'double');

% Generate X-axis data as a sequence of integers from 1
x = 1:numel(data);

% Plot the entanglement entropy
plot(x, data, 'o'); hold on;
xlabel('Site','Interpreter','latex')
ylabel('Entanglement Entropy','Interpreter','latex')
set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);


% Display the variable values used in the filename
disp(['L = ', num2str(L)]);
disp(['Delta = ', num2str(Delta)]);
disp(['h = ', num2str(h)]);
disp(['gamma = ', num2str(gamma)]);
disp(['D = ', num2str(D)]);

if dump_figure
    % Save the figure as EPS file with the specified name
    figure_name_eps = ['ee_brownianL', num2str(L), 'Delta', num2str(Delta), ...
        'h', num2str(h), 'gamma', num2str(gamma), 'D', num2str(D), '.eps'];
    figure_path = fullfile('../note_figure', figure_name_eps);
    saveas(gcf, figure_path, 'epsc');
    disp(['Figure saved as: ', figure_path]);
end

%{
figure;
fitting_formular = 'FMH'; % 'FMH' for Ferromagnetic Heisenberg ground state; 'CFT'

start_site = L/4;
if(strcmp(fitting_formular, 'FMH') )
    plot(log( x .*(L-x)/L),data,'-o');hold on;
    x_select = x(start_site: 1: end-start_site);
    fit_x=log( x_select.*(L-x_select)/L);
    xlabel_str = '$\log(x (L-x)/L)$';
elseif(strcmp(fitting_formular, 'CFT'))
    plot(log(sin(pi*(x/L))),data,'-o');hold on;
    fit_x=log(sin(pi.*(x(start_site:1:end-start_site)/L)));
    xlabel_str = '$\log(\sin(\pi x /L))$';
end

fit_y=data(start_site:1:end-start_site);
p = fit((fit_x'),(fit_y),'poly1');
plot(fit_x, p.p1*fit_x + p.p2,'-.');hold on;
c = p.p1 * 6;
fprintf(' c =%.5f\n',c);

xlabel(xlabel_str,'Interpreter','latex')
ylabel('Entanglement Entropy','Interpreter','latex')
set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);
%}