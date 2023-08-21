%support z and xyz channel
lambda = 0:0.04:0.5;
L_list = [60; 64; 68; 124;128;132; 252;256;260];
channel_type = 'z';
log_g_data = zeros(numel(L_list)/3, numel(lambda));
L_central_list = zeros(numel(L_list)/3, 1);
for i = 1: numel(L_list)/3
    L_small = L_list(3*i-2); L_central = L_list(3*i-2); L_large = L_list(3*i);
    L_central_list(i) = L_central;
    for j = 1:numel(lambda)
        file_name = ['../data/renyi2_entropyL', num2str(L_small), 'channel', channel_type, 'lambda',num2str(lambda(j), '%.6f')];
        file_id = fopen(file_name,'r');
        S2_small = fread(file_id,1, 'double');
        fclose(file_id);

        file_name = ['../data/renyi2_entropyL', num2str(L_central), 'channel', channel_type, 'lambda',num2str(lambda(j), '%.6f')];
        file_id = fopen(file_name,'r');
        S2 = fread(file_id,1, 'double');
        fclose(file_id);

        file_name = ['../data/renyi2_entropyL', num2str(L_large), 'channel', channel_type, 'lambda',num2str(lambda(j), '%.6f')];
        file_id = fopen(file_name,'r');
        S2_large = fread(file_id,1, 'double');
        fclose(file_id);

        log_g_data(i,j) = -(S2 - (S2_large - S2_small)/8 * L_central);
    end
end

plot(lambda, log_g_data,'-o');hold on;

% plot(L_list', thermal_renyi2','-o');hold on;

l=legend('$L=64$', '$128$', '$256$');
set(l,'Box','off');set(l,'Interpreter','latex');
set(l,'Fontsize',24);
set(l,'Location','SouthWest');

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$\lambda$','Interpreter','latex');
ylabel('$\log(g)$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);

figure;
plot(log(L_central_list), (log_g_data(:,end)),'-o');hold on;
set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds

figure;

for i = 1: numel(L_central_list)
    L = L_central_list(i);
    plot(lambda*L^(0.3), log_g_data(i,:),'-o');hold on;
end


set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$\lambda$','Interpreter','latex');
ylabel('$\log(g)$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);