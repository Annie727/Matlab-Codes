function [deltaVals, gammaVals, muVals] = PEParamsFromData_1 ...
            (inputCaseInfo, ProvinceName, dateCuts)

% Getting specific columns of the case information file
opts = detectImportOptions(inputCaseInfo);
opts.SelectedVariableNames = {'Province', ...
    'Date_of_Report', ...
    'Date_of_Onset', ...
    'Case_Classification', ...
    'Outcome', ...
    'Date_Recovered', ...
    'Date_Died', ...
    'Vaccination_Status'};

% Reading the file
T = readtable(inputCaseInfo, opts);

% Capitalizing string columns to make the data uniform
T.Province = upper(T.Province);
T.Case_Classification = upper(T.Case_Classification);
T.Outcome = upper(T.Outcome);
T.Vaccination_Status = upper(T.Vaccination_Status);

% We want to find the proportion of V that becomes E (delta),
% recovery rate (gamma) and death rate (mu)
deltaVals = zeros(1, length(dateCuts) - 1);
gammaVals = zeros(1, length(dateCuts) - 1);
muVals = zeros(1, length(dateCuts) - 1);

for i = 1:length(dateCuts) - 1
    startCut = datetime(dateCuts{i}, 'Format', 'MM/dd/yyyy');
    endCut = datetime(dateCuts{i + 1}, 'Format', 'MM/dd/yyyy');

    if i == 1
        CQidx = find(T.Date_of_Report <= endCut);
    else
        CQidx = find(T.Date_of_Report > startCut & T.Date_of_Report <= endCut);
    end

    % proportion of V that becomes E (delta)
    Vidx = find(contains(T(CQidx, :).Vaccination_Status, 'FULLY VACCINATED'));
    deltaList = length(Vidx);
%     disp(deltaList);
    if deltaList == 0
        delta = 0;
        deltaVals(i) = delta;
    else
        kappa=1084;
%         disp(mean(deltaList));
        delta = deltaList/kappa;
        deltaVals(i) = delta;
    end
    
    % recovery rate (gamma)
    gammaList = days(T(CQidx, :).Date_Recovered ...
                            - T(CQidx, :).Date_of_Onset);
    gammaList = gammaList(gammaList >= 0);
    gamma = mean(gammaList);
    gammaVals(i) = gamma;

    % death rate
    if all(contains(T(CQidx, :).Date_Died, 'N/A'))
        mu = Inf;
        muVals(i) = mu;
    else
        dateDied = datetime(T(CQidx, :).Date_Died, ...
                    "Format", "dd-MMM-yyyy", "InputFormat", "dd-MMM-yyyy");
        muList = days(dateDied - T(CQidx, :).Date_of_Onset);
        muList = muList(muList >= 0);
        mu = mean(muList);
        muVals(i) = mu;
    end
end

gammaVals = max(0, min(1, 1 ./ gammaVals));
muVals = max(0, min(1, 1 ./ muVals));
end