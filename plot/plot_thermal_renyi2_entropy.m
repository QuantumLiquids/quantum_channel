% each row corresponds to different v
% each column corresponds to different L
L_list = 10:10:80;
v_list = 0:0.1:0.5;
g_list = zeros(1, numel(v_list));
select_L = 1:1:numel(L_list);

thermal_renyi2_entropy = [  0  0   0    0   0  0      0         0; ...
    1.7451    3.4809    5.2129    6.9429    8.6717   10.3997   12.1271   13.8540;...
    2.8426    5.6530    8.4507   11.2418   14.0289   16.8132   19.5954   22.3762;...
    3.3183    6.5873    9.8369   13.0767   16.3103   19.5399   22.7664   25.9907;...
    3.4346    6.8148   10.1736   13.5216   16.8630   20.1998   23.5334   26.8645;...
    3.4424    6.8301   10.1963   13.5516   16.9002   20.2443   23.5850   26.9233];
% (thermal_renyi2_entropy(:,1:end-1) + thermal_renyi2_entropy(:,2:end))/2- diff(thermal_renyi2_entropy, 1);
for i = 1:numel(v_list)
    v = v_list(i);
    plot(L_list, thermal_renyi2_entropy(i,:),'-o');hold on;
    p = fit((L_list(select_L)'),thermal_renyi2_entropy(i,select_L)','poly1');
    fprintf('v = %.2f, log(g)=%.5f\n',v, p.p2);
    g_list(i) = exp(p.p2);
end


set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$L$','Interpreter','latex');
ylabel('Thermal renyi-2 entropy $S$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);

l=legend('$\lambda=0$', '$0.1$','$0.2$','$0.3$','$0.4$','$0.5$');
set(l,'Box','off');set(l,'Interpreter','latex');
set(l,'Fontsize',24);
set(l,'Location','SouthWest');


figure;
plot(v_list, log(g_list),'-o'); hold on;
set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$\lambda$','Interpreter','latex');
ylabel('$\log(g)$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);

    