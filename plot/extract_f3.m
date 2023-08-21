% each row corresponds to different v
% each column corresponds to different L
L_list = 10:10:80;
v_list = 0:0.1:0.5;
f3_list = zeros(1, numel(v_list));
select_L = 2: 2:numel(L_list);
thermal_renyi3_entropy = [   0       0         0    0         0         0         0         0;...     
    1.4504    2.8867    4.3177    5.7455    7.1714    8.5957   10.0189   11.4411;...
    2.5360    4.9550    7.3343    9.6928   12.0384   14.3752   16.7056   19.0311;...
    3.0191    5.8041    8.5288   11.2254   13.9056   16.5752   19.2371   21.8935;...
    3.1302    5.9917    8.7900   11.5593   14.3120   17.0537   19.7877   22.5160;...
    3.1374    6.0037    8.8068   11.5807   14.3380   17.0843   19.8229   22.5559
];
negativity3 = [0.7035    0.2693    0.7340    0.3507    0.7558    0.3994    0.7726    0.4342;...
    0.3843    0.1873    0.3743    0.2302    0.3705    0.2520    0.3682    0.2658;...
    0.1234    0.0636    0.0923    0.0654    0.0815    0.0651    0.0761    0.0646;...
    0.0204    0.0101    0.0132    0.0100    0.0117    0.0099    0.0110    0.0098;...
    0.0011    0.0005    0.0007    0.0005    0.0006    0.0005    0.0006    0.0005;...
         0         0         0         0         0         0         0         0];


for i = 1:numel(v_list)
    v = v_list(i);
    plot(log(L_list(select_L)), negativity3(i,(select_L)),'-o');hold on;
    p = fit(log(L_list(select_L)'),negativity3(i,(select_L))','poly1');
    fprintf('lambda = %.2f, f2=%.5f\n',v, p.p1);
    f3_list(i) = p.p1;
end


set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$log(L)$','Interpreter','latex');
ylabel('$S_A^{(3)}$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);


figure;
plot(v_list, f3_list,'-o'); hold on;
set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$\lambda$','Interpreter','latex');
ylabel('$f_3$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);
