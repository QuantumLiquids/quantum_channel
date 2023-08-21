% raw data only keep 5th digital
lambda = 0:0.02:0.5;
L_list = [60; 68; 92;100; 140; 148;188;196;508;516];
load("thermal_renyi2_Heisenberg_add_Zchannel",'S2')
thermal_renyi2=S2;
log_g_data = zeros(numel(L_list)/2, numel(lambda));
for i = 1: numel(L_list)/2
    L_small = L_list(2*i-1); L_large = L_list(2*i);
    Delta_L = L_large - L_small;
    L_central = (L_small + L_large)/2;
    log_g_data(i,:) = (thermal_renyi2(2*i,:) + thermal_renyi2(2*i-1,:))/2 - (thermal_renyi2(2*i,:) - thermal_renyi2(2*i-1,:))/Delta_L * L_central;
end
log_g_data = -log_g_data;
plot(lambda, log_g_data,'-o');hold on;

% plot(L_list', thermal_renyi2','-o');hold on;

l=legend('$L=64$', '$96$','$144$','$196$','$512$');
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
L_list_central = [64, 96, 144, 196, 512];
semilogx((L_list_central), (log_g_data(:,6)),'-o');hold on;
set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds

figure;

for i = 1: numel(L_list_central)
   L = L_list_central(i);
   plot(lambda*L^(0.3), log_g_data(i,:),'-o');hold on;
end


set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$\lambda$','Interpreter','latex');
ylabel('$\log(g)$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);