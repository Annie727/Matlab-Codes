function dataXLSXtoMAT_1(inputCaseInfo, ProvinceName, ...
    inputProvinceData, dateCuts, days2test)

% Getting specific columns of the case information file
opts = detectImportOptions('Region_XII.xlsx');
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

% Comparing start and end of date cuts with the first and last days of
% report, respectively
uniqueReport = unique(T.Date_of_Report);
startDate = datetime(dateCuts{1}, 'Format', 'MM/dd/yyyy');
endDate = datetime(dateCuts{end}, 'Format', 'MM/dd/yyyy') + days2test;

if startDate ~= uniqueReport(1)
    error(append('Error. Change the start date of dateCuts to ', string(uniqueReport(1)), ...
        ' to match with the first date of report.'))
end
if endDate ~= uniqueReport(end)
    error(append('Error. Change the end date of dateCuts to ', string(uniqueReport(end - 31)), ...
        ' to match with the last date of report.'))
end

% Tabulating the time-series data
dates = startDate:endDate;
dataDate = dates;
dataCase = zeros(1, length(dates)).';
dataDeath = zeros(1, length(dates)).';
dataRecov = zeros(1, length(dates)).';

for i = 1:length(dates)
    dateidx = find(T.Date_of_Report == dates(i));

    dataCase(i) = length(find(contains(T(dateidx, :).Case_Classification, 'CONFIRMED')));

    dataDeath(i) = length(find(contains(T(dateidx, :).Outcome, 'DIED')));

    dataRecov(i) = length(find(contains(T(dateidx, :).Outcome, 'RECOVERED')));
end

% And saving the file
save(inputProvinceData, 'dataDate', 'dataCase', 'dataDeath', 'dataRecov');

end