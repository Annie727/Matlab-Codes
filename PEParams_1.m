function PEParams_1(inputCaseInfo, ProvinceName, dateCuts, ...
            inputProvinceParams)

% Fixed parameter values
beta = 0.1; % transmission rate, placeholder
sigma = 50; % multiplicative, placeholder
alpha = 0.2; % transition rate or incubation period
rho = 0.1; % reporting rate
phi = 1/7; % severity rate
tau = 0.1; % proportion of Iu that becomes Ir, placeholder
epsilon = 0.05; % vaccine inefficacy rate
kappa = 1084; % average number of fully vaccinated individuals per day
lambda = 1/180; % vaccine waning rate
theta = 1/240; % immunity loss rate
zeta = 1/26017.2; % natural death rate
% Zeta is from https://population.un.org/wpp/Download/Standard/Population/

% Getting parameter values based from data
currDir = cd;
% cd('data');
[deltaVals, gammaVals, muVals] = PEParamsFromData_1(inputCaseInfo, ...
            ProvinceName, dateCuts);

pset = zeros(14, length(deltaVals));

% Assign parameters to vector element
% Note that the parameter set is a 3D matrix since we have parameters for
% each CQ/date cut
pset(1, :) = beta .* ones(1, length(dateCuts) - 1);
pset(2, :) = sigma .* ones(1, length(dateCuts) - 1);
pset(3, :) = alpha .* ones(1, length(dateCuts) - 1);
pset(4, :) = rho .* ones(1, length(dateCuts) - 1);
pset(5, :) = phi .* ones(1, length(dateCuts) - 1);
pset(6, :) = tau .* ones(1, length(dateCuts) - 1);
pset(7, :) = epsilon .* ones(1, length(dateCuts) - 1);
pset(8, :) = kappa .* ones(1, length(dateCuts) - 1);
pset(9, :) = lambda .* ones(1, length(dateCuts) - 1);
pset(10, :) = deltaVals;
pset(11, :) = gammaVals;
pset(12, :) = theta .* ones(1, length(dateCuts) - 1);
pset(13, :) = muVals;
pset(14, :) = zeta .* ones(1, length(dateCuts) - 1);

% And saving the parameter file
save(inputProvinceParams, 'pset');
cd(currDir)

end