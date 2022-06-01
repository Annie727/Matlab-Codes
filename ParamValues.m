global beta sigma alpha rho phi tau epsilon kappa lambda delta gamma theta mu zeta;
ICs = [975476, 1, 0, 0, 0, 0, 0, 0, 0];
tspan = 0:1:275;

beta = 0.0136;
sigma = 15.6313;
alpha = 0.2000;
rho = 0.1000;
phi = 0.1429;
tau = 0.0634;
epsilon = 0.4000;
kappa = 5000;
lambda = 1/180;
delta = 0.1503;
gamma = 0.0481;
theta = 1/240;
mu = 0.0685;
zeta = 1/26017.2;

[time, solution] = ode15s(@CovidModel, tspan, ICs);

figure(1)
plot(tspan, solution(:,3)', 'DisplayName', 'confirmed cases')
hold on;

plot(tspan, solution(:,5)', 'DisplayName', 'recoveries')
hold on;

plot(tspan, solution(:,8)', 'DisplayName', 'deaths')

lgd = legend;
lgd.FontSize = 14;

xlim([0, 275]);
xlabel('Time in Days');
ylabel('Count');

title('Simulation 3')
txt = '{\kappa} = 5000, {\epsilon} = 0.01';
subtitle(txt)

saveas(gcf, 'Simulation3.fig')
exportgraphics(gcf, 'Simulation3.png', 'Resolution', 300)