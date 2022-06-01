function dstate = ODEModel_1(~,state,param,a)

% Assign values to the parameters
beta = param(1);
sigma = param(2); 
alpha = param(3);
rho = param(4);
phi = param(5);
tau = param(6);
epsilon = param(7);
kappa = param(8);
lambda = param(9);
delta = param(10);
gamma = param(11);
theta = param(12);
mu = param(13);

zeta = a(1);

% beta = param(1);
% sigma = param(2); 
% tau = param(3);
% delta = param(4);
% gamma = param(5);
% mu = param(6);
% 
% alpha = a(1);
% rho = a(2);
% phi = a(3);
% epsilon = a(4);
% kappa = a(5);
% lambda = a(6);
% theta = a(7);
% zeta = a(8);

% Assign values to the variables
S = state(1);
E = state(2);
Ir = state(3);
Iu = state(4);
R = state(5);
V = state(6);
D = state(7);
N = S + E + Ir + Iu + R + V;

% Model Equations
dSdt = theta*R ...
       +lambda*(1-delta)*V ...
       -beta*S*(Ir+sigma*Iu)/N ...
       -kappa ...
       -zeta*S;
dEdt = (S+epsilon*delta*V)*beta*(Ir+sigma*Iu)/N ...
       -(alpha+zeta)*E;
dIrdt = alpha*rho*E ...
        +phi*tau*Iu ...
        -(gamma+mu+zeta)*Ir;
dIudt = alpha*(1-rho)*E ...
        -(phi*tau+gamma*(1-tau)+zeta)*Iu;
dRdt = gamma*(Ir+(1-tau)*Iu) ...
       -(theta+zeta)*R;
dVdt = kappa ...
       -lambda*(1-delta)*V ...
       -beta*epsilon*delta*V*(Ir+sigma*Iu)/N ...
       -zeta*V;
dDdt = mu*Ir;

% Create vector for state variables
dstate = zeros(7,1);

% Assign the equations to the model variables
dstate(1) = dSdt;
dstate(2) = dEdt;
dstate(3) = dIrdt;
dstate(4) = dIudt;
dstate(5) = dRdt;
dstate(6) = dVdt;
dstate(7) = dDdt;

end            