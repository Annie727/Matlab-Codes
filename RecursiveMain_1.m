% CODE FOR PARAMETER ESTIMATION, BOOTSTRAPPING, AND FORECASTING

clear
clc
close all

LHSsamples = 10;
BSsamples = 10;

%Output summary detailing the run, an appendix to naming LHS and BS files
outputSummary = ['Final5'];

% Case information XLSX file inside \data
inputCaseInfo = 'Region_XII.xlsx';

% Parameters to be estimated
param2est = {'beta', 'sigma', 'tau'};
% Data to be used (default is cases)
data2use = {'C'}; 
% Number of days to be tested (meaning, last day of training is last day of
% report minus days2test
days2test = 7;
% Number of days to forecast (this already includes days2test)
days2FC = 31;

FuncSouthCotabato(LHSsamples, BSsamples, outputSummary, inputCaseInfo, param2est, data2use, days2test, days2FC);