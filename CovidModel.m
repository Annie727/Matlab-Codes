function dstate = CovidModel(~,state)

global beta sigma alpha rho phi tau epsilon kappa lambda delta gamma theta mu zeta;

% assign values to the parameters
S = state(1);
E = state(2);
Ir = state(3);
Iu = state(4);
R = state(5);
V = state(6);
C = state(7);
D = state(8);
Rec = state(9);
N = S+E+Ir+Iu+R+V;

% model equations
dSdt = 975476/26017.2 ...
       +theta*R ...
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
dCdt = alpha*rho*E+phi*tau*Iu;
% dDdt = (alpha*rho*E+phi*tau*Iu)*mu;
dDdt = zeta*N+mu*Ir;
dRecdt = (alpha*rho*E+phi*tau*Iu)*gamma;

% create vector of length 9x1
dstate = zeros(9,1);

% assign elements to the model which are the differential equations in the
% model
dstate(1) = dSdt;
dstate(2) = dEdt;
dstate(3) = dIrdt;
dstate(4) = dIudt;
dstate(5) = dRdt;
dstate(6) = dVdt;
dstate(7) = dCdt;
dstate(8) = dDdt;
dstate(9) = dRecdt;

end