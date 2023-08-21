% L = 128, xyz channel, Renyi-2 mutual information
% I^(2)(A,B) = S^(2)_A + S^(2)_B - S^(2)_{AB},  S^(2)_{AB} is the thermal
% renyi entropy

L = 256;
channel_type = 'z';
lambda_set = 0: 0.04:0.48;
La_set = 1:L-1;
x = log(sin(La_set/L * pi));
select_data = 2+L/4:2:L-1-L/4;
f2_list = zeros(1, numel(lambda_set));
for i = 1:numel(lambda_set)
    lambda = lambda_set(i);
    file_name = ['../data/renyi2_entropyL', num2str(L), 'channel', channel_type, 'lambda',num2str(lambda, '%.6f')];
    file_id = fopen(file_name,'r');
    thermal_entropy = fread(file_id,1, 'double');
    entangle_a = fread(file_id,L-1, 'double');
    entangle_b = fread(file_id,L-1, 'double');
    fclose(file_id);
    mutual_information_data = entangle_a + flip(entangle_b) - thermal_entropy;  % column vector
    % plot(x(select_data), mutual_information_data(select_data),'-o'); hold on;

    p = fit(x(select_data)',mutual_information_data(select_data),'poly1');
    fprintf('lambda = %.2f, f_2=%.5f\n',lambda, p.p1/2);
    f2_list(i) = p.p1/2;
end

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$\log(\sin(\pi L_A/L))$','Interpreter','latex');
ylabel('$I^{(2)}(A,\bar A)$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);

% figure;

plot(lambda_set, f2_list, '-o');hold on;

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$\lambda$','Interpreter','latex');
ylabel('$f_2$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);