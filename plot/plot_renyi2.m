L = 256;
file_name = ['../data/renyi2_entropyL',num2str(L), 'channelzlambda0.000000'];
file_id = fopen(file_name,'r');
thermal_entropy = fread(file_id,1, 'double'); %0.0
entangle_a = fread(file_id,L-1, 'double');
entangle_b = fread(file_id,L-1, 'double');
ee2 = entangle_a;

la_list = 1:L-1;

start_site = 4;

plot(la_list, ee2,'-o');hold on;


modelfun = @(b,x)(b(1)/8 * log(sin(pi*(2 .* x + 1)/2./(L+1))) + b(2) - b(3) * sin(pi/2*(2.*x+1))./sqrt(sin(pi*(2 * x + 1)/2./(L+1))) );
mdl = fitnlm(la_list(start_site:1:end-start_site+1),ee2(start_site:1:end-start_site+1),modelfun,[1,0.7,sqrt(pi/L)]);

c = mdl.Coefficients.Estimate(1);
En = mdl.Coefficients.Estimate(2);
% kF = mdl.Coefficients.Estimate(4);
b = mdl.Coefficients.Estimate;
sites = la_list(start_site):0.1:la_list(end-start_site+1);
plot(sites, modelfun(b,sites),'-');


% plot(log((L+1)/pi * sin(pi*(2*l_list+1)/2/(L+1))),ee2,'-o');hold on;
% fit_x=log((L+1)/pi * sin(pi*(2*l_list(start_site:2:end-start_site)+1)/2/(L+1)));
% fit_y=ee2(start_site:2:end-start_site);
% p = fit((fit_x'),(fit_y'),'poly1');

% plot(fit_x, p.p1*fit_x + p.p2,'-.');hold on;
fprintf('c =%.5f\n',c);

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$\log(\sin(l/L))$','Interpreter','latex');
ylabel('Renyi 2 EE','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);
