function FuncSouthCotabato(LHSsamples, BSsamples, outputSummary, inputCaseInfo, param2est, data2use, days2test, days2FC)
tic;
disp('STARTING SOUTH COTABATO')
ProvinceName = 'SOUTH COTABATO';

% Population of the area being considered
areaPop = 975476;

% Date cuts for each CQ or arbitrary cut to be considered
dateCuts = {'03/01/2021', ... % START
            '04/01/2021', ... % Arbitrary cut
            '05/01/2021', ... % Arbitrary cut
            '06/01/2021', ... % Arbitrary cut
            '06/13/2021', ... % Arbitrary cut
            '07/01/2021', ... % Arbitrary cut
            '08/01/2021', ... % Arbitrary cut
            '09/01/2021', ... % Arbitrary cut
            '09/16/2021', ... % Arbitrary cut
            '10/01/2021', ... % Arbitrary cut
            '11/01/2021', ... % Arbitrary cut
            '11/30/2021'}; % Arbitrary cut

% Name of each date cut
dateCQStr = {'', '', '', '', '', '', '', '', '', '', ''};

% Did you make changes to above conditions?
% true - reruns data file and parameter file
changes2Conds = true;

% Names of data and parameter files based on case information file name
inputProvinceData = append('Data_',ProvinceName,'_',...
                    erase(inputCaseInfo, '.xlsx'), '.mat');
inputProvinceParams = append('Params_',ProvinceName,'_',...
                    erase(inputCaseInfo, '.xlsx'), '.mat');

% Processing data file
% If-condition checks if you have already processed the data file so that
% you don't have to run it again to save time
currDir = cd;
% cd('data');
if ~isfile(inputProvinceData) || changes2Conds == true
    disp('Processing data file from XLSX file.')
    dataXLSXtoMAT_1(inputCaseInfo, ProvinceName, inputProvinceData, ...
        dateCuts, days2test);
else
    disp('Data file from XLSX file already present.')
end
disp('Data file done.')
disp('=====================');


% Processing parameter file
% Same with data file, if-condition checks if the parameter file is already
% present
if ~isfile(inputProvinceParams) || changes2Conds == true
    disp('Processing parameter file from XLSX file.')
    cd(currDir)
    PEParams_1(inputCaseInfo, ProvinceName, dateCuts, ...
            inputProvinceParams);
else
    disp('Parameter file from XLSX file already present.')
end
disp('Parameter file done.')
disp('=====================');
cd(currDir)

% Parameter estimation with Latin hypercube sampling
LHSNameResults = ParameterEstimation_1(inputProvinceData, inputProvinceParams, ProvinceName, ...
    areaPop, dateCuts, outputSummary, param2est, data2use, LHSsamples, days2test);

% Parameter bootstrapping with short-term forecasting
Bootstrapping_1(LHSNameResults, inputProvinceData, inputProvinceParams, ProvinceName, ...
    areaPop, dateCuts, outputSummary, param2est, data2use, BSsamples, ...
    dateCQStr, days2FC, days2test);
end