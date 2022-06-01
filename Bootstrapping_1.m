function Bootstrapping_1(LHSNameResults, inputProvinceData, inputProvinceParams, ProvinceName, ...
    areaPop, dateCuts, outputSummary, param2est, data2use, BSsamples, ...
    dateCQStr, days2FC, days2test)
tic;

inputLHSResults = append(LHSNameResults, '.mat');

maData = true;
maDays = 10;

popFactor = 1;

E0 = 100;

%% i) Training data

currDir = cd; % Access data in the 'data' folder
% cd('data');
dataDiv = load(inputProvinceData);
paramsetFromFile = load(inputProvinceParams);
paramsetParams = paramsetFromFile.pset;
cd(currDir)

cd('lhs_results_1');
LHSFullResults = load(inputLHSResults);
LHSoutputCurves = LHSFullResults.outputCurves;
LHSoutputIncid = LHSFullResults.outputIncid;
LHSoutputPSet = LHSFullResults.outputPSet;
LHSoutputP1 = LHSFullResults.outputP1;
LHSoutputRes = LHSFullResults.outputRes;
LHSoutputMet = LHSFullResults.outputMet;
LHSoutputModSelect = LHSFullResults.outputModSelect;
LHSoutputRepNo = LHSFullResults.outputRepNo;
LHSoutputResid = LHSFullResults.outputResid;
cd(currDir);

SSEVals = sum(sum(LHSoutputRes(:, 1:2, :), 3), 2);
% disp('SSEVals');
% disp(SSEVals);
[~, minG] = min(SSEVals);
p0Vals = LHSoutputP1(minG, :, :);
p0Vals = reshape(p0Vals, [length(dateCuts) - 1, length(param2est)]);


% Start and end of preferred whole period
startDate = find(dataDiv.dataDate == datetime(dateCuts{1}, 'Format', 'MM/dd/yyyy'));
endDate = find(dataDiv.dataDate == datetime(dateCuts{end}, 'Format', 'MM/dd/yyyy'));
endDateFC = find(dataDiv.dataDate == datetime(dateCuts{end}, 'Format', 'MM/dd/yyyy') + days2test);

dataDate = dataDiv.dataDate;
dataCase = double(dataDiv.dataCase(startDate:endDate));
dataDeath = double(dataDiv.dataDeath(startDate:endDate));
dataRecov = double(dataDiv.dataRecov(startDate:endDate));

dataCaseMA = movmean(dataCase, maDays);
dataDeathMA = movmean(dataDeath, maDays);
dataRecovMA = movmean(dataRecov, maDays);

dataCaseFC = double(dataDiv.dataCase(startDate:endDateFC));
dataCaseMAFC = movmean(dataCaseFC, maDays);

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
Iu0 = 100; %* (areaPop/975476);
R0 = data2EstR(1);
V0 = 100;                                                     
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
LHSmatrix = Model_LHS_1(LB, UB, BSsamples, 'unif', 1e20); % LHS call function
cd(currDir)

lenDateCuts = length(dateCuts) - 1;

% Pre-allocation of cells/matrices
inputData = [dataCase dataDeath dataRecov]; 
inputEst = [data2EstC data2EstD data2EstR]; 
inputP0 = LHSmatrix;

outputIncid = cell(BSsamples, lenDateCuts);
outputCurves = cell(BSsamples, lenDateCuts);
outputData = cell(BSsamples, lenDateCuts);
outputPSet = zeros(BSsamples, 22, lenDateCuts);
outputP1 = zeros(BSsamples, length(paramIdx), lenDateCuts);
outputRes = zeros(BSsamples, 1 * 4, lenDateCuts);
outputMet = zeros(BSsamples, 4, lenDateCuts);
outputRepNo = zeros(BSsamples, lenDateCuts);
outputResid = cell(BSsamples, length(data2use), lenDateCuts);
outputModSelect = zeros(BSsamples, 4, lenDateCuts);

outputFCInits = zeros(BSsamples, 8);

popReset = pop;

disp('Start of parameter bootstrapping.')

for g = 1:BSsamples
    disp(append('   BS: ', string(g), ' out of ', string(BSsamples)));
    for j = 1:lenDateCuts
        p0 = p0Vals(j, :);
%         disp(p0);

        options = optimoptions('lsqcurvefit', ...
                                'Algorithm', 'trust-region-reflective', ...
                                'MaxFunctionEvaluations', 1e4*length(p0), ...
                                'MaxIterations', 1e4*length(p0), ...
                                'Display', 'none');

        paramset(1:14) = paramsetParams(:, j);

        % Daily data
        dataDivC = LHSoutputIncid{minG, j}(:, 1);

        data4Est = zeros(length(tPeriods{j}), length(data2use));
        for i = 1:length(data2use)
            incid2Poiss = eval(append('dataDiv', data2use{i}));
            caseData = zeros(length(incid2Poiss), 1);
            caseData(1) = incid2Poiss(1);
            for tt = 2:length(incid2Poiss)
                caseData(tt, 1) = poissrnd(incid2Poiss(tt), 1, 1);
            end 
            data4Est(:, i) = caseData;
        end
%         disp('data4Est');
%         disp(data4Est);
        outputData{g, j} = data4Est;
%         disp('outputData');
%         disp(outputData);

        xDataLCF = repmat(tPeriods{j}, 1, length(data2use));
        yDataLCF = reshape(data4Est, 1, numel(data4Est));
%         disp('xDataLCF');
%         disp(xDataLCF);
%         disp('yDataLCF');
%         disp(yDataLCF);
%         disp('Start of PEObjectiveLCF');
        [P, fval, ~, exitflag, output] = lsqcurvefit(@PEObjectiveLCF_1, ...
            p0, xDataLCF, yDataLCF, LB, UB, options, ...
            paramIdx, paramset, pop, data4Est, data2use);
%         disp('End of PEObjectiveLCF');
        
        % For-loop placing the estimated values of parameters
        for i = 1:length(paramIdx)
            paramset(paramIdx(i)) = P(i);
        end
        outputP1(g, :, j) = P;
        outputPSet(g, :, j) = paramset;
%         disp('outputP1');
%         disp(outputP1);
%         disp('outputPSet');
%         disp(outputPSet);

        initCmpts = paramset(15:22);
%         disp('initCmpts');
%         disp(initCmpts);
        % Solving the ODE again using the parameter estimates and adjustments
%         disp('Start of BaselineModel');
        odeOptions = odeset('Reltol', 1e-6, 'Abstol', 1e-6);
        [~, sol] = ode15s(@BaselineModel_1, tPeriods{j}, initCmpts, ...
            odeOptions, paramset);
        outputCurves{g, j} = sol;
%         disp('End of BaselineModel');

        Csol = sol(:, 7);

        Cdiff = [Csol(1); diff(Csol)].';
        incidAll = [Cdiff.'];
        outputIncid{g, j} = incidAll;

        for i = 1:8
            paramset(i + 14) = sol(end, i);
        end
        paramset(21) = Cdiff(end);
%         disp('paramset(15:22)');
%         disp(paramset(15:22));

        pop = sol(end, 1) + sol(end, 2) + sol(end, 3) ...
                + sol(end, 4) + sol(end, 5) + sol(end, 6);

%         disp('pop')
%         disp(pop);

        outputFCInits(g, :) = paramset(15:22);

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

        modSelectArray = PEModelSelect_1(dataDivC, ...
                         param2est, SSEc);
        outputModSelect(g, :, j) = modSelectArray;
    end
    paramset = paramsetReset;
    pop = popReset;
end
disp('End of parameter bootstrapping.')
toc;

outputFCIncid = cell(BSsamples, lenDateCuts + 3);
outputFCCurves = cell(BSsamples, lenDateCuts + 3);
outputFCPSet = zeros(BSsamples, 22, lenDateCuts + 3);
outputFCRepNo = zeros(BSsamples, lenDateCuts + 3);

outputFCIncid(:, 1:lenDateCuts) = outputIncid;
outputFCCurves(:, 1:lenDateCuts) = outputCurves;
outputFCPSet(:, :, 1:lenDateCuts) = outputPSet;
outputFCRepNo(:, 1:lenDateCuts) = outputRepNo;

for i = 1:lenDateCuts
    beta(i) = plims_1(outputP1(:,1 , i), 0.5);
end

betaMedian = median(beta, 'all');
betaAve = mean(beta, 'all');
betaLast = outputP1(:, 1, lenDateCuts);
betaMid = plims_1(betaLast, 0.5);

beta4FC = [betaMedian betaMid betaAve];

j = lenDateCuts;
tFCPeriod = (tPeriods{j}(end) + 1):(tPeriods{j}(end) + days2FC);
for k = 1:3
%     disp(k);
    for g = 1:BSsamples
        paramset = outputPSet(g, :, j);
        paramset(1) = beta4FC(k);
        paramset(15:22) = outputFCInits(g, :);
        outputFCPSet(g, :, j + k) = paramset;
    
        initCmpts = paramset(15:22);
        % Solving the ODE again using the parameter estimates and adjustments
        odeOptions = odeset('Reltol', 1e-6, 'Abstol', 1e-6);
        [~, sol] = ode15s(@BaselineModel_1, tFCPeriod, initCmpts, ...
            odeOptions, paramset);
        outputFCCurves{g, j + k} = sol;
    
        Csol = sol(:, 7);
    
        Cdiff = [Csol(1); diff(Csol)].';
        
        incidAll = [Cdiff.'];
        outputFCIncid{g, j + k} = incidAll;
    
        for i = 1:8
            paramset(i + 13) = sol(end, i);
        end
        paramset(21) = Cdiff(end);
    
        RepNo = PERepNo_1(paramset);
        outputFCRepNo(g, j + k) = RepNo;

    end
end


% Exporting data determined
BSNameResults = append('Result_BS_',ProvinceName,'_',string(BSsamples),'_',...
                            string(today('datetime')), '_', ...
                             outputSummary);

disp('File name of results:');
disp(BSNameResults);
disp('=====================');

simsDir = 'bs_results_1';
cd(simsDir);
save(append(BSNameResults, '.mat'), ...
    'tPeriods', 'timeVect', 'inputData', 'inputEst', 'inputP0', ...
    'outputCurves', 'outputData', 'outputIncid', 'outputPSet', ...
    'outputP1', 'outputRes', 'outputMet', ...
    'outputModSelect', 'outputRepNo', 'outputResid', ...
    'outputFCIncid', 'outputFCCurves', 'outputFCPSet', ...
    'outputFCRepNo','beta', 'betaLast','betaMid',...
    'betaMedian', 'betaAve');
cd(currDir);

BSPlots_1(tPeriods, timeVect, dataDate, inputData, inputEst, ...
    outputCurves, outputData, outputIncid, outputP1, outputRepNo, ...
    areaPop, BSsamples, dateCQStr, days2test, BSNameResults, ProvinceName, outputSummary)

FCPlots_1(tPeriods, timeVect, dataDate, inputData, inputEst, ...
    minG, areaPop, BSsamples, days2FC, outputFCIncid, outputFCCurves, ...
    outputFCRepNo, beta4FC, days2test, dataCaseFC, ...
    dataCaseMAFC, BSNameResults, ProvinceName, outputSummary)

end