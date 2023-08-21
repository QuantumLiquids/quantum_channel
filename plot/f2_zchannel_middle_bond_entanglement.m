% each row corresponds to different v
% each column corresponds to different L
L_list = 10:10:80;
v_list = 0:0.1:0.5;
f2_list = zeros(1, numel(v_list));
select_L = 2:2:numel(L_list);
thermal_renyi2_entropy = [  0  0   0    0   0  0      0         0; ...
    1.7451    3.4809    5.2129    6.9429    8.6717   10.3997   12.1271   13.8540;...
    2.8426    5.6530    8.4507   11.2418   14.0289   16.8132   19.5954   22.3762;...
    3.3183    6.5873    9.8369   13.0767   16.3103   19.5399   22.7664   25.9907;...
    3.4346    6.8148   10.1736   13.5216   16.8630   20.1998   23.5334   26.8645;...
    3.4424    6.8301   10.1963   13.5516   16.9002   20.2443   23.5850   26.9233];
ee2 = [  0.7068    0.3470    0.7472    0.4467    0.7756    0.5052    0.7973    0.5467;...
    1.4024    2.0311    3.1681    3.8423    4.9202    5.6172    6.6652    7.3769;...
    1.8413    3.0789    4.6783    5.9444    7.4905    8.7714   10.2912   11.5820;...
    2.0322    3.5298    5.3281    6.8434    8.5899   10.1154   11.8368   13.3694;...
    2.0790    3.6396    5.4863    7.0617    8.8568   10.4412   12.2113   13.8021;...
    2.0821    3.6470    5.4970    7.0765    8.8748   10.4631   12.2366   13.8312];

log_term = ee2 - thermal_renyi2_entropy/2;


for i = 1:numel(v_list)
    v = v_list(i);
    plot(log(L_list(select_L)), log_term(i,(select_L)),'-o');hold on;
    p = fit(log(L_list(select_L)'),log_term(i,(select_L))','poly1');
    fprintf('v = %.2f, f2=%.5f\n',v, p.p1);
    f2_list(i) = p.p1;
end


set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$log(L)$','Interpreter','latex');
ylabel('$S_A^{(2)} - S^{(2)}/2$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);


figure;
plot(v_list, f2_list,'-o'); hold on;
set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$v$','Interpreter','latex');
ylabel('$f_2$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);
