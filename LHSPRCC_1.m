%Susan Massey
%Most General version of the LHS code

clear
close all
clc

%%_________________________________________________________________________
%
% FIRST STEP: Latin Hypercube Sampling of Model Parameters 
%__________________________________________________________________________

% This requests user input via the command line to get started:
M = input('Number of parameters to sample?: '); % Don't include add'l parameters you would like to leave fixed
while rem(M,1)~=0 || M<=0
    M = input('Number of parameters should be an integer > 0. Please re-enter the number of parameters: ');
end

N = input('Number of samples to draw? (recommend 100 to 1000): '); 
while rem(N,1)~=0 || N<=0
    N = input('Number of samples should be an integer > 0. Please re-enter the number of samples: ');
end

% This code will prompt for sample distribution specifics and return drawn and randomly paired parameter samples
[parameters,A] = DrawSamples_1(M,N);

% EDIT THE FOLLOWING FOR SAVING YOUR SAMPLES AND SIMULATION DATA:
outFileStr = 'LHS-ODE'; % give workspace an appropriate unique name
outFileName1 = [outFileStr,'_samples.mat'];
outFileName2 = [outFileStr,'_results.mat'];


save(outFileName1, 'parameters', 'A')

%---Plot Histograms to Visualize Parameter Distributions---% 
 
% Specify how many bins you want for your histograms:                            
binCount = ceil(N/100); % can change this manually, just make sure it's an integer value  

histPlot = plotSampleHists_1(M,parameters,binCount);

pause(3) %Time to dock/maximize the figure before it saves, if you prefer.

% Name and save histograms visualizing the sample distributions:
figurelabel1=([outFileStr,'-N',num2str(N),'-histograms.fig']);
figurelabel2=([outFileStr,'-N',num2str(N),'-histograms.png']);
saveas(histPlot, figurelabel1);
saveas(histPlot, figurelabel2);
%%
%_________________________________________________________________________
%
% SECOND STEP: Solve Model Function with Sampled Parameters
%__________________________________________________________________________

% EDIT THE FOLLOWING VARIABLES, UNSAMPLED PARAMETERS, & ANY OTHER ARGS HERE

% Specify independent var(s) to pass to the model function
x = 0:1:275; % time from March 8, 2020 to March 5, 2021

unsampledps(1) = 1/26017.2; %zeta

% South Cotabato province parameters on March 2021
pop = 975476;
E0 = 100;
Ir0 = 1;
Iu0 = 0;
R0 = 0;
V0 = 0;
D0 = 0;
S0 = pop-E0-Iu0-Ir0-R0-V0;

Init=[S0 E0 Iu0 Ir0 R0 V0 D0];

PRCCSimdata = struct; % initialize struct to store outputs for computing PRCCs

tic % start measuring time to solve equations to monitor progress

% Loop over the parameter sample pairs for Monte Carlo simulations
for j=1:N   
    sampledparams=A(j,1:M);
    fprintf('Parameters passed for sample: ');fprintf('%u',j);fprintf(' of ');fprintf('%u\n',N);
    
    odeoptions = odeset('Reltol',1e-12,'Abstol',1e-12);
    [t,sol] = ode15s(@ODEModel_1,x,Init,odeoptions,sampledparams,unsampledps);
    PRCCSimdata(j).c = sol(:,3); % outcomemeasure 1
    PRCCSimdata(j).u = sol(:,4); % outcomemeasure 2
    
    % Assign parameters to compute R0
    beta = sampledparams(1);
    sigma = sampledparams(2); 
    alpha = sampledparams(3);
    rho = sampledparams(4);
    phi = sampledparams(5);
    tau = sampledparams(6);
    epsilon = sampledparams(7);
    kappa = sampledparams(8);
    lambda = sampledparams(9);
    delta = sampledparams(10);
    gamma = sampledparams(11);
    theta = sampledparams(12);
    mu = sampledparams(13);

    zeta = unsampledps(1);
    
    suscept = sol(:,1);

    % Basic Rep. No equation:
    W = alpha+zeta; %A
    X = gamma+mu+zeta;  %B
    Y = gamma+zeta-gamma*tau+phi*tau;  %C

    RepNo = (alpha*beta*(gamma*rho+phi*tau+rho*zeta-gamma*rho*tau ...
        -sigma*(rho-1)*X))/(W*X*Y);
    Rt = (RepNo.*suscept)./pop;
    
    PRCCSimdata(j).r = Rt; % outcomemeasure 3
    toc
end

save(outFileName2, 'PRCCSimdata') % add Simdata to the saved .mat file

OutputOfInterest_Ir = zeros(N,length(x));
OutputOfInterest_Iu = zeros(N,length(x));
OutputOfInterest_Rt = zeros(N,length(x));

for si = 1:N
% EDIT THE FOLLOWING TO SPECIFY OUTPUT DATA TO COMPARE:

    OutputOfInterest_Ir(si,:) = PRCCSimdata(si).c;
    OutputOfInterest_Iu(si,:) = PRCCSimdata(si).u;
    OutputOfInterest_Rt(si,:) = PRCCSimdata(si).r;
end


% EDIT THE FOLLOWING STRING TO NAME OUTPUT DATA:
labelstring_1 = 'I_r cases'; % consistent id for plot labels, filenames
labelstring_2 = 'I_u cases';
labelstring_3 = 'R_t';

%---Visualize the range of simulation results with errorbar plots---%

outPlot_1 = plotSimulationOutput_1(x,N,OutputOfInterest_Ir,labelstring_1); 
outPlot_2 = plotSimulationOutput_1(x,N,OutputOfInterest_Iu,labelstring_2); 
outPlot_3 = plotSimulationOutput_1(x,N,OutputOfInterest_Rt,labelstring_3); 

% SUMMARIZE
OutputOfInterest = {OutputOfInterest_Ir,OutputOfInterest_Iu,OutputOfInterest_Rt};
labelstring = {labelstring_1,labelstring_2,labelstring_3};
outPlot = {outPlot_1,outPlot_2,outPlot_3};

pause(7)

for i = 1:length(OutputOfInterest)
    figurelabel1=([outFileStr,'-N',num2str(N),'-',labelstring{i},'-ErrorbarPlot.fig']);
    figurelabel2=([outFileStr,'-N',num2str(N),'-',labelstring{i},'-ErrorbarPlot.png']);
    saveas(outPlot{i},figurelabel1);
    saveas(outPlot{i},figurelabel2);
end

%%
%%_________________________________________________________________________
%
% THIRD STEP: Compute & Plot Partial Rank Correlation Coefficients (PRCC)
%__________________________________________________________________________
prcc_col = cell(1,length(OutputOfInterest));
prccPlot = cell(1,length(OutputOfInterest));
studentT_col = cell(1,length(OutputOfInterest));

for i = 1:length(OutputOfInterest)
% EDIT INDICES FOR EVALUATING RESULTS BELOW:

    % To look at variation in SPACE:

    t_index =[]; % for a spatially varied ODE model
    [prcc,studentT] = VariedPRCC_1(M,N,A,OutputOfInterest{i},x,t_index);
    prcc_col{i} = prcc;
    studentT_col{i} = studentT;
    prccPlot{i} = plotVariedPRCC_1(M,N,x,labelstring{i},parameters,prcc_col{i});

end

save('PRCCvalues.mat', 'prcc_col') 
save('studentTstats.mat', 'studentT_col') 

pause(5) %Time to dock/maximize the figure before it saves, if you prefer.
for i = 1:length(OutputOfInterest)
    figurelabel1=([outFileStr,'-N',num2str(N),'-', labelstring{i},'-PRCC.fig']); %assign an appropriate filename
    figurelabel2=([outFileStr,'-N',num2str(N),'-', labelstring{i},'-PRCC.png']);
    saveas(prccPlot{i},figurelabel1); %save fig
    saveas(prccPlot{i},figurelabel2); %save png
end

%%
% Extract prcc for specific time points
p=categorical({parameters.name});
prccPlot_spec = cell(1,length(OutputOfInterest));
col_time = cell(1,length(OutputOfInterest));

for i = 1:length(OutputOfInterest)
    prccvar_col = prcc_col{1,i};
    k1 = prccvar_col(:,1);
    k2 = prccvar_col(:,162);
    k3 = prccvar_col(:,218);
    k4 = prccvar_col(:,276);
    
    prcctime = [k1, k2, k3,k4];
    col_time{i} = prcctime;
    
    figure()
    hold on
    box on
    bar(p,prcctime(1:13,:));
    newcolors = [0.50 0.65 0.15;0 0.5 1; 0.5 0 1; 0.7 0.7 0.7];
    colororder(newcolors)
    
    ylabel('PRCC value');
    legend({'End of Q1','End of Q2','End of Q3','End of Q4'},'Location','EastOutside','Fontsize',15)
    title(['Correlation Plot of ',labelstring{i},' from LHS simulations, ',num2str(N),' samples']);
    prccPlot_spec{i} = gca;
end
pause(5) %Time to dock/maximize the figure before it saves, if you prefer.

for i = 1:length(OutputOfInterest)
	figurelabel1=([outFileStr,'-N',num2str(N),'-', labelstring{i},'-PRCCunvar.fig']); %assign an appropriate filename
    figurelabel2=([outFileStr,'-N',num2str(N),'-', labelstring{i},'-PRCCunvar.png']);
    saveas(prccPlot_spec{i},figurelabel1); %save fig
    saveas(prccPlot_spec{i},figurelabel2); %save png
end