function [dstate] = BaselineModel_1(t, state, param)
% disp('Start of Baseline Model');
% COVID-19 model with vaccination
% disp('Start of BaselineModel')
% Assign values to the parameters
beta = param(1);
% disp('beta')
% disp(beta);
sigma = param(2); 
alpha = param(3);
% disp('alpha');
% disp(alpha);
rho = param(4);
phi = param(5);
tau = param(6);
epsilon = param(7);
% disp('epsilon');
% disp(epsilon);
kappa = param(8);
lambda = param(9);
delta = param(10);
gamma = param(11);
theta = param(12);
mu = param(13);
zeta = param(14);
% disp('zeta');
% disp(zeta);

% Assign values to the variables
% disp('Assign values to the variables');
S = state(1);
% disp('S');
% disp(S);
E = state(2);
% disp('E');
% disp(E);
Ir = state(3);
% disp('Ir');
% disp(Ir);
Iu = state(4);
% disp('Iu');
% disp(Iu);
R = state(5);
% disp('R');
% disp(R);
V = state(6);
% disp('V');
% disp(V);
C = state(7);
% disp('C');
% disp(C);
D = state(8);
% disp('D');
% disp(D);
N = S + E + Ir + Iu + R + V;
% disp('N');
% disp(N);

% Model Equations
% disp('Model Equations');
dSdt = theta*R ...
       +lambda*(1-delta)*V ...
       -beta*S*(Ir+sigma*Iu)/N ...
       -kappa ...
       -zeta*S;
% disp('dSdt');
% disp(dSdt);
dEdt = (S+epsilon*delta*V)*beta*(Ir+sigma*Iu)/N ...
       -(alpha+zeta)*E;
% disp('dEdt');
% disp(dEdt);
dIrdt = alpha*rho*E ...
        +phi*tau*Iu ...
        -(gamma+mu+zeta)*Ir;
% disp('dIrdt');
% disp(dIrdt);
dIudt = alpha*(1-rho)*E ...
        -(phi*tau+gamma*(1-tau)+zeta)*Iu;
% disp('dIudt');
% disp(dIudt);
dRdt = gamma*(Ir+(1-tau)*Iu) ...
       -(theta+zeta)*R;
% disp('R');
% disp(R);
dVdt = kappa ...
       -lambda*(1-delta)*V ...
       -beta*epsilon*delta*V*(Ir+sigma*Iu)/N ...
       -zeta*V;
% disp('V');
% disp(V);
dCdt = alpha*rho*E ...
       +phi*tau*Iu;
% disp('C');
% disp(C);
dDdt = mu*Ir;
% disp('dDdt');
% disp(dDdt);


% Create vector for state variables
dstate = zeros(8,1);

% Assign the equations to the model variables
dstate(1) = dSdt;
dstate(2) = dEdt;
dstate(3) = dIrdt;
dstate(4) = dIudt;
dstate(5) = dRdt;
dstate(6) = dVdt;
dstate(7) = dCdt;
dstate(8) = dDdt;

% disp(dstate(1));
% disp('End of BaselineModel');