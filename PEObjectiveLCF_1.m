% A Function that will solve the ODE model to obtain cumulative cases/incidence

function APPENDSOLN = PEObjectiveLCF_1(p0,xDataLCF,idxparams,paramset,pop,casesEst,data2use)
% disp('Start of PEObjectiveLCF');
% disp('p0');
% disp(p0);
% disp('xDataLCF');
% disp(xDataLCF);
% disp('length(p0)');
% disp(length(p0));
% Assign the initial values to the parameters that needed to be estimated
k=1;

for i=1:length(p0)
    paramset(idxparams(k)) = p0(i);
%     disp(k);
%     disp('paramset(idxparams(k))');
%     disp(paramset(idxparams(k)));
%     disp('p0(i)');
%     disp(i);
%     disp(p0(i));
    k=k+1;
end
% disp('End of for loop');
% Assign the values (from the parameter set) to variable names here
E0 = paramset(16);
% disp('E0');
% disp(E0);
Ir0 = paramset(17);
% disp('Ir0');
% disp(Ir0);
Iu0 = paramset(18);
% disp('Iu0');
% disp(Iu0);
R0 = paramset(19);
% disp('R0');
% disp(R0);
V0 = paramset(20);
% disp('V0');
% disp(V0);
C0 = paramset(21);
% disp('C0');
% disp(C0);
D0 = paramset(22);
% disp('D0');
% disp(D0);

S0 = pop - (E0 + Ir0 + Iu0 + R0 + V0 + D0);
% disp('S0');
% disp(S0);

init = [S0 E0 Ir0 Iu0 R0 V0 C0 D0];
% disp('init');
% disp(init);

% timeVect = xDataLCF(1:(length(xDataLCF)/length(data2use)));
timeVect = xDataLCF(1:(length(xDataLCF)));
% disp('timeVect');
% disp(timeVect);

options = odeset('Reltol',1e-12,'Abstol',1e-12);
% disp(options);
% disp('ode45');
[~, sol] = ode45(@BaselineModel_1,timeVect,init,options,paramset);

% Defining the output of the least-squares function
APPENDSOLN = zeros(size(casesEst));
modelSoln = sol(:, 7);
modelDiff = [modelSoln(1); diff(modelSoln)].';
APPENDSOLN(:, 1) = modelDiff;

APPENDSOLN = reshape(APPENDSOLN, 1, numel(APPENDSOLN));
% disp(APPENDSOLN);
% disp('End of PEObjectiveLCF');
end