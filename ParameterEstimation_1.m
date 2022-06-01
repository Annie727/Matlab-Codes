function LHSNameResults = ParameterEstimation_1(inputProvinceData, inputProvinceParams, ProvinceName, ...
    areaPop, dateCuts, outputSummary, param2est, data2use, LHSsamples, days2test)
tic;

maData = true;
maDays = 10;

popFactor = 1;

E0 = 100;
Iu0 = 0; %* (areaPop/975476);
V0 = 0;

%% i) Training data

currDir = cd; % Access data in the 'data' folder
% cd('data');
dataDiv = load(inputProvinceData);
paramsetFromFile = load(inputProvinceParams);
paramsetParams = paramsetFromFile.pset;
cd(currDir)

% Start and end of preferred whole period
startDate = find(dataDiv.dataDate == datetime(dateCuts{1}, 'Format', 'MM/dd/yyyy'));
endDate = find(dataDiv.dataDate == datetime(dateCuts{end}, 'Format', 'MM/dd/yyyy'));
% disp(startDate);
% disp(endDate);
dataDate = dataDiv.dataDate;
dataCase = double(dataDiv.dataCase(startDate:endDate));
dataDeath = double(dataDiv.dataDeath(startDate:endDate));
dataRecov = double(dataDiv.dataRecov(startDate:endDate));

dataCaseMA = movmean(dataCase, maDays);
dataDeathMA = movmean(dataDeath, maDays);
dataRecovMA = movmean(dataRecov, maDays);
% disp(dataCaseMA);
tCuts = [0];
tPeriods = cell(1, length(dateCuts) - 1);
for i = 1:length(dateCuts)-1
    tCuts(end + 1) = days(datetime(dateCuts{i+1}, 'InputFormat', 'MM/dd/yyyy')...
                        - datetime(dateCuts{1}, 'InputFormat', 'MM/dd/yyyy'));
    tPeriods{i} = tCuts(i):1:tCuts(i+1);
end
timeVect = startDate-1:1:endDate-1;


%% Parameters and their Bounds
 
% parameter names in the order they are displayed in the results file
paramNames = {'beta', 'sigma', 'alpha', 'rho', 'phi', ...
                'tau', 'epsilon', 'kappa', 'lambda', ...
                'delta', 'gamma', 'theta', 'mu', 'zeta', ...
                'S0', 'E0', 'Ir0', 'Iu0', 'R0', 'V0', 'C0', 'D0'};

% Pre-allocate for the LB&UB of param2est
paramIdx = zeros(1,length(param2est));

pop = areaPop*popFactor;
if maData == true
    data2EstC = dataCaseMA;
    data2EstD = dataDeathMA;
    data2EstR = dataRecovMA;
else
    data2EstC = dataCase;
    data2EstD = dataDeath;
    data2EstR = dataRecov;
end

Ir0 = data2EstC(1);                                                                       
R0 = data2EstR(1);                                                      
C0 = data2EstC(1);
D0 = data2EstD(1);
S0 = pop - E0 - Ir0 - Iu0 - R0 - V0 - D0;

paramset = zeros(22, 1);
paramset(1:14) = paramsetParams(:, 1);
for i = 15:22
    paramset(i) = eval(paramNames{i});
end
paramsetReset = paramset;

% bounds for the parameters and variables
LBfull = [1e-5,1,1e-5,1e-5,1e-5,1e-5,1e-5,1,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5];
UBfull = [1,100,1,1,1,1,1,1e7,1,1,1,1,1,1,1e7,1e7,1e7,1e7,1e7,1e7,1e7,1e7];
    
LB = zeros(1, length(param2est));
UB = zeros(1, length(param2est));

% Assign the LB and UB to the parameters in the param2est
for i = 1:length(param2est)
    k = 1;
    while strcmp(param2est(i),paramNames(k)) == 0 %e.g. false/not the same param name
        k = k + 1;
    end
    paramIdx(i) = k;
    LB(i) = LBfull(k);
    UB(i) = UBfull(k);
end

dependDir = 'dependencies';
cd(dependDir)
LHSmatrix = Model_LHS_1(LB, UB, LHSsamples, 'unif', 1e20); % LHS call function
cd(currDir)

lenDateCuts = length(dateCuts) - 1;

% Pre-allocation of cells/matrices
inputData = [dataCase dataDeath dataRecov]; 
inputEst = [data2EstC data2EstD data2EstR]; 
inputP0 = LHSmatrix;

outputIncid = cell(LHSsamples, lenDateCuts);
outputCurves = cell(LHSsamples, lenDateCuts);
outputPSet = zeros(LHSsamples, 22, lenDateCuts);
outputP1 = zeros(LHSsamples, length(paramIdx), lenDateCuts);
outputRes = zeros(LHSsamples, 1 * 4, lenDateCuts);
outputMet = zeros(LHSsamples, 4, lenDateCuts);
outputRepNo = zeros(LHSsamples, lenDateCuts);
outputResid = cell(LHSsamples, length(data2use), lenDateCuts);
outputModSelect = zeros(LHSsamples, 4, lenDateCuts);

popReset = pop;

SSELHS = zeros(1, LHSsamples);

disp('Start of parameter estimation.')

for g = 1:LHSsamples
    p0 = LHSmatrix(g,:);

    options = optimoptions('lsqcurvefit', ...
                            'Algorithm', 'trust-region-reflective', ...
                            'MaxFunctionEvaluations', 1e4*length(p0), ...
                            'MaxIterations', 1e4*length(p0), ...
                            'FunctionTolerance', 1e-10, ...
                            'OptimalityTolerance', 1e-10, ...
                            'StepTolerance', 1e-10, ...
                            'Display', 'none');
    SSEg = 0;
    disp(append('   LHS: ', string(g), ' out of ', string(LHSsamples)));
    for j = 1:lenDateCuts
        paramset(1:14) = paramsetParams(:, j);
        k = 1;
        for l = 1:length(p0)
            paramset(paramIdx(k)) = p0(l);
            k = k+1;
        end

        % Daily data
        dataDivC = data2EstC(tPeriods{j}+1);

        data4Est = zeros(length(tPeriods{j}), length(data2use));
        for i = 1:length(data2use)
            data4Est(:, i) = eval(append('dataDiv', data2use{i}));
        end
%         disp('data4Est');
%         disp(data4Est);
%         data4Est = zeros(length(tPeriods{j}), 1);
%         data4Est(:, 1) = eval(append('dataDiv', data2use{1}));
%         disp('data4Est');
%         disp(data4Est);

        xDataLCF = tPeriods{j};
        yDataLCF = reshape(data4Est, 1, numel(data4Est));
%         disp('p0');
%         disp(p0);
%         disp('xDataLCF');
%         disp(xDataLCF);
%         disp('yDataLCF');
%         disp(yDataLCF);
%         display(paramset)
        [P, fval, ~, exitflag, output] = lsqcurvefit(@PEObjectiveLCF_1, ...
            p0, xDataLCF, yDataLCF, LB, UB, options, ...
            paramIdx, paramset, pop, data4Est, data2use);
        
        % For-loop placing the estimated values of parameters
        for i = 1:length(paramIdx)
            paramset(paramIdx(i)) = P(i);
        end
        outputP1(g, :, j) = P;
        outputPSet(g, :, j) = paramset;

        initCmpts = paramset(15:22);
        % Solving the ODE again using the parameter estimates and adjustments
        odeOptions = odeset('Reltol', 1e-6, 'Abstol', 1e-6);
        [~, sol] = ode15s(@BaselineModel_1, tPeriods{j}, initCmpts, ...
            odeOptions, paramset);
        outputCurves{g, j} = sol;

        Csol = sol(:, 7);

        Cdiff = [Csol(1); diff(Csol)].';
        incidAll = [Cdiff.'];
        outputIncid{g, j} = incidAll;

        for i = 1:8
            paramset(i + 14) = sol(end, i);
        end

        paramset(21) = Cdiff(end);

        pop = pop - sol(end, 8);

        RepNo = PERepNo_1(paramset);
        outputRepNo(g, j) = RepNo;


        % Residuals
        residC = Cdiff - dataDivC.';
        outputResid{g, 1, j} = residC;

        SSEc  = sum((Cdiff - dataDivC.').^2);
        RMSEc = sqrt(mean((Cdiff - dataDivC.').^2)); %root-mean-squared error
        MAEc  = mean(abs(Cdiff - dataDivC.')); %mean absolute error
        MAPEc = mean(abs((Cdiff - dataDivC.')./dataDivC.')); %mean absolute percentage error

        resultsArray = [SSEc,...
                        RMSEc, ...
                        MAEc, ...
                        MAPEc];
        metricsArray = [exitflag, fval, ...
                        output.funcCount, output.iterations];
        outputRes(g, :, j) = resultsArray;
        outputMet(g, :, j) = metricsArray;

        SSEg = SSEg + SSEc;

        modSelectArray = PEModelSelect_1(dataDivC, ...
                         param2est, SSEc);
        outputModSelect(g, :, j) = modSelectArray;
    end
    SSELHS(g) = SSEg;
    paramset = paramsetReset;
    pop = popReset;
end
disp('End of parameter estimation.')
toc;

[~, minG] = min(SSELHS);

% Exporting data determined
LHSNameResults = append('Result_LHS_',ProvinceName,'_',string(LHSsamples),'_',...
                            string(today('datetime')), '_', ...
                             outputSummary);

disp('File name of results:');
disp(LHSNameResults);
disp('=====================');

simsDir = 'lhs_results_1';
cd(simsDir);
save(append(LHSNameResults, '.mat'), ...
    'tPeriods', 'timeVect', 'inputData', 'inputEst', 'inputP0', ...
    'outputCurves', 'outputIncid', 'outputPSet', 'outputP1', 'outputRes', 'outputMet', ...
    'outputModSelect', 'outputRepNo', 'outputResid');
cd(currDir);

PEPlots_1(tPeriods, timeVect, dataDate, inputData, inputEst, ...
    outputCurves, outputIncid, outputRepNo, minG, areaPop, LHSNameResults, ...
    days2test)

end