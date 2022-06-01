function [ModSelect] = PEModelSelect_1(dataDivC, param2est, SSEc)

% Akaike Information Criteria
Nc = length(dataDivC); %No. of data points

s = length(param2est); %No. of parameters fitted 

% AIC Rule of Thumb: n/K >= 40
Rule_c = Nc/s; %For confirmed cases

% AIC equation
AIC_c = Nc*log(SSEc/Nc)+2*(s+1); %For confirmed cases

% AICc equation e.g. if n/K < 40
AICc_c = AIC_c + (2*s*(s+1))/(Nc-s-1); %For confirmed cases

% Bayesian Information Criteria
BIC_c = Nc*log(SSEc/Nc)+ (s+1)*log(Nc)/Nc; 

Rule = [Rule_c];
AIC = [AIC_c];
AICc = [AICc_c];
BIC = [BIC_c];

% Store results to a new cell
ModSelect = [Rule, AIC, AICc, BIC];

end