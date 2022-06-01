function BSPlots_1(tPeriods, timeVect, dataDate, inputData, inputEst, ...
    outputCurves, outputData, outputIncid, outputP1, ...
    outputRepNo, cityPop, BSsamples, dateCQStr, days2test, BSNameResults, ProvinceName, outputSummary)

currDir = cd;
figsDir = 'bs_figs_1';
cd(figsDir)
if ~exist(BSNameResults, 'dir')
   mkdir(BSNameResults)
end
cd(BSNameResults)

color_model = 1/255*[30,144,255];
color_cases = 1/255*[255,160,122];
color_casesra = 1/255*[255,127,80];
color_S = 1/255*[255,218,185];
color_E = 1/255*[218,165,32];
color_C = 1/255*[255,160,122];
color_Iu = 1/255*[255,238,130];
color_R = 1/255*[144,238,144];
color_V = 1/255*[218,165,32];                                      % !! ??
color_death = 1/255*[112,128,144];
color_boot = 1/255*[75,0,130];

startMo = dateshift(dataDate(1), 'start', 'month');
endMo = dateshift(dataDate(end-days2test), 'start', 'month');
if startMo < dataDate(1)
    startMo = startMo + calmonths(1);
    duraMo = split(between(startMo, endMo, 'months'), 'months');
    dateXLabels = sort([startMo + calmonths(0:1:duraMo)]);
    dateX = datefind(dateXLabels, dataDate);
    xDate = datestr(dateXLabels, 'mmm yyyy');
else
    duraMo = split(between(startMo, endMo, 'months'), 'months');
    dateXLabels = sort([startMo + calmonths(0:1:duraMo)]);
    dateX = datefind(dateXLabels, dataDate);
    xDate = datestr(dateXLabels, 'mmm yyyy');
end

CDAll = zeros(BSsamples, length(timeVect));
CAll = zeros(BSsamples, length(timeVect));
SAll = zeros(BSsamples, length(timeVect));
EAll = zeros(BSsamples, length(timeVect));
IrAll = zeros(BSsamples, length(timeVect));
IuAll = zeros(BSsamples, length(timeVect));
RAll = zeros(BSsamples, length(timeVect));
VAll = zeros(BSsamples, length(timeVect));
DAll = zeros(BSsamples, length(timeVect));
R0All = zeros(BSsamples, length(timeVect));
for g = 1:BSsamples
    for i = 1:length(tPeriods)
        CAll(g, tPeriods{i} + 1) = outputIncid{g, i}(:, 1);
        CDAll(g, tPeriods{i} + 1) = outputData{g, i}(:, 1);

        SAll(g, tPeriods{i} + 1) = outputCurves{g, i}(:, 1);
        EAll(g, tPeriods{i} + 1) = outputCurves{g, i}(:, 2);
        IrAll(g, tPeriods{i} + 1) = outputCurves{g, i}(:, 3);
        IuAll(g, tPeriods{i} + 1) = outputCurves{g, i}(:, 4);
        RAll(g, tPeriods{i} + 1) = outputCurves{g, i}(:, 5);
        VAll(g, tPeriods{i} + 1) = outputCurves{g, i}(:, 6);
        DAll(g, tPeriods{i} + 1) = outputCurves{g, i}(:, 8);

        Ssol = outputCurves{g, i}(:, 1);
        Nsol = cityPop - outputCurves{g, i}(:, 8);
        R0All(g, tPeriods{i} + 1) = outputRepNo(g, i) .* (Ssol ./ Nsol);
    end
end


% Daily data vs fit
figure(5)
tiledlayout(1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
bar(timeVect, inputData(:, 1), 'FaceColor', color_cases, ...
    'EdgeColor', 'none')
hold on
plot(timeVect, inputEst(:, 1), 'color', color_casesra, 'LineWidth', 4);
hold on

cd(currDir)
CMed = plims_1(CAll, 0.5);
CUp = plims_1(CAll, 0.975);
CLow = plims_1(CAll, 0.025);
cd(figsDir)
cd(BSNameResults)

casesFit = plot(timeVect, CMed, '-', 'Color', color_model, ...
    'LineWidth', 2);
hold on

fillX = [timeVect, fliplr(timeVect)];
fillY = [CLow, fliplr(CUp)];
h = fill(fillX, fillY, color_model, 'linestyle', 'none');
set(h,'facealpha', .5)
hold on

tPStart = zeros(1, length(tPeriods));
for i = 1:length(tPeriods)
    xline(tPeriods{i}(1), 'color', color_death, 'LineStyle', '--', ...
        'LineWidth', 2)
    tPStart(i) = tPeriods{i}(1);
    hold on
end

ylabel('new confirmed cases', 'FontSize', 18)
set(gca, 'XTick', tPStart, ...
         'XTickLabel', dateCQStr, ...
         'XTickLabelRotation', 45)
set(gca, 'FontSize', 15)

xlabel('dates', 'FontSize', 18)
set(gca, 'XTick', dateX - 1, ...
         'XTickLabel', xDate, ...
         'XTickLabelRotation', 45)
set(gca, 'FontSize', 15)

saveas(gcf, 'BSDaily.fig')
exportgraphics(gcf, 'BSDaily.png', 'Resolution', 300)

% Cumulative data vs fit
figure(6)
tiledlayout(1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
bar(timeVect, cumsum(inputData(:, 1)), 'FaceColor', color_cases, ...
    'EdgeColor', 'none')
hold on
plot(timeVect, cumsum(inputEst(:, 1)), 'color', color_casesra, 'LineWidth', 4);
hold on

cd(currDir)
CMed = plims_1(cumsum(CAll, 2), 0.5);
CUp = plims_1(cumsum(CAll, 2), 0.975);
CLow = plims_1(cumsum(CAll, 2), 0.025);
cd(figsDir)
cd(BSNameResults)

casesFit = plot(timeVect, CMed, '-', 'Color', color_model, ...
    'LineWidth', 2);
hold on

fillX = [timeVect, fliplr(timeVect)];
fillY = [CLow, fliplr(CUp)];
h = fill(fillX, fillY, color_model, 'linestyle', 'none');
set(h,'facealpha', .5)
hold on

tPStart = zeros(1, length(tPeriods));
for i = 1:length(tPeriods)
    xline(tPeriods{i}(1), 'color', color_death, 'LineStyle', '--', ...
        'LineWidth', 2)
    tPStart(i) = tPeriods{i}(1);
    hold on
end

ylabel('cumulative confirmed cases', 'FontSize', 18)
set(gca, 'XTick', tPStart, ...
         'XTickLabel', dateCQStr, ...
         'XTickLabelRotation', 45)
set(gca, 'FontSize', 15)

saveas(gcf, 'BSCumulative.fig')
exportgraphics(gcf, 'BSCumulative.png', 'Resolution', 300)

% Compartments
figure(7)
tiledlayout(6, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
cmpts = {'S', 'E', 'C', 'Iu', 'R', 'V'};
cmpts1 = {'S', 'E', 'Ir', 'Iu', 'R', 'V'};
compartments = {'susceptible', 'exposed', 'reported', ...
            'unreported', 'recovered', 'vaccinated'};

for j = 1:length(cmpts)
%     disp('cmpts')
%     disp(cmpts{j});
    nexttile;

    cd(currDir)
    XMed = plims_1(eval(append(cmpts1{j}, 'All')), 0.5);
    XUp = plims_1(eval(append(cmpts1{j}, 'All')), 0.975);
    XLow = plims_1(eval(append(cmpts1{j}, 'All')), 0.025);
    cd(figsDir)
    cd(BSNameResults)

    casesFit = plot(timeVect, XMed, '-', 'Color', eval(append('color_', cmpts{j})), ...
        'LineWidth', 2);
    hold on
    
    fillX = [timeVect, fliplr(timeVect)];
    fillY = [XLow, fliplr(XUp)];
%     if j==4
%         disp(fillX);
%         disp(fillY);
%     end

    h = fill(fillX, fillY, eval(append('color_', cmpts{j})), 'linestyle', 'none');
    set(h,'facealpha', .5)
    hold on
    
    for i = 1:length(tPeriods)
        xline(tPeriods{i}(1), 'color', color_death, 'LineStyle', '--', ...
            'LineWidth', 2)
        hold on
    end

    ylabel(compartments{j}, 'FontSize', 18)
    xlim([0 tPeriods{i}(end)])
    set(gca, 'FontSize', 15)
    set(gca, 'XTick', [])
end
xlabel('dates', 'FontSize', 18)
set(gca, 'XTick', dateX - 1, ...
         'XTickLabel', xDate, ...
         'XTickLabelRotation', 45)
set(gca, 'FontSize', 15)

saveas(gcf, 'BSCompartments.fig')
exportgraphics(gcf, 'BSCompartments.png', 'Resolution', 300)

% Time-varying reproduction number
figure(8)
tiledlayout(1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;

cd(currDir)
R0Med = plims_1(R0All, 0.5);
R0Up = plims_1(R0All, 0.975);
R0Low = plims_1(R0All, 0.025);
cd(figsDir)
cd(BSNameResults)

casesFit = semilogy(timeVect, R0Med, '-', 'Color', color_model, ...
    'LineWidth', 2);
grid on
hold on

R0BSCI(:, 1) = [R0Med];
R0BSCI(:, 2) = [R0Low];
R0BSCI(:, 3) = [R0Up];

fillX = [timeVect, fliplr(timeVect)];
fillY = [R0Low, fliplr(R0Up)];
h = fill(fillX, fillY, color_boot, 'linestyle', 'none');
set(h,'facealpha', .5)
hold on

for i = 1:length(tPeriods)
    xline(tPeriods{i}(1), 'color', color_death, 'LineStyle', '--', ...
        'LineWidth', 2)
    hold on
end

ylabel('time-varying reproduction number', 'FontSize', 18)
xlim([0 tPeriods{i}(end)])

yline(1, 'color', color_death, ...
        'LineStyle', '-', 'LineWidth', 4);

xlabel('dates', 'FontSize', 18)
set(gca, 'XTick', dateX - 1, ...
         'XTickLabel', xDate, ...
         'XTickLabelRotation', 45)
set(gca, 'FontSize', 15)

saveas(gcf, 'BSRepNo.fig')
exportgraphics(gcf, 'BSRepNo.png', 'Resolution', 300)


% Histogram of estimates
figure(9)
tiledlayout(3, length(tPeriods), 'Padding', 'compact', 'TileSpacing', 'compact');

for i = 1:3
    for j = 1:length(tPeriods)
        nexttile
        PE2H = outputP1(:, i, j);
    
        cd(currDir)
        PEMed = plims_1(PE2H, 0.5);
        PEUp = plims_1(PE2H, 0.975);
        PELow = plims_1(PE2H, 0.025);
% disp(i)
% disp(j)        
% display(PEMed)
% display(PEUp)
% display(PELow)
    if i==1
        betaCI(:, 3*j - 2) = [PEMed];
        betaCI(:, 3*j - 1) = [PELow];
        betaCI(:, 3*j) = [PEUp];
    elseif i==2
        sigmaCI(:, 3*j - 2) = [PEMed];
        sigmaCI(:, 3*j - 1) = [PELow];
        sigmaCI(:, 3*j) = [PEUp];
    else
        tauCI(:, 3*j - 2) = [PEMed];
        tauCI(:, 3*j - 1) = [PELow];
        tauCI(:, 3*j) = [PEUp];
    end
        cd(figsDir)
        cd(BSNameResults)

        xline(PEMed, '-', 'Color', color_boot, ...
            'LineWidth', 2, ...
            'Label', {[num2str(round(PEMed, 4))];...
                        ['[', num2str(round(PELow, 4))];...
                        ['-', num2str(round(PEUp, 4)), ']']}, ...
            'LabelOrientation', 'horizontal');
        hold on
        
        fillX = [PELow, PELow, PEUp, PEUp];
        fillY = [0, 1000, 1000, 0];

        h = fill(fillX, fillY, color_boot, 'linestyle', 'none');
        set(h,'facealpha', 0.2)
        hold on

        histplot = histogram(PE2H, 'EdgeColor', 'white', ...
            'FaceColor', color_boot, 'NumBins', 10);
        ylim([0, max(histplot.Values) + mod(-max(histplot.Values), BSsamples)])
    
        if j == 1 && i == 1
            ylabel('\beta');
        elseif j == 1 && i == 2
            ylabel('\sigma');
        elseif j==1 && i==3
            ylabel('\tau');
        end
        set(gca, 'FontSize', 15)
    end
end

saveas(gcf, 'BSHistEsts.fig')
exportgraphics(gcf, 'BSHistEsts.png', 'Resolution', 300)

BS2NameResults = append('Result_BS2_',ProvinceName,'_',string(BSsamples),'_',...
                            string(today('datetime')), '_', ...
                             outputSummary);

disp('File name of results:');
disp(BS2NameResults);
disp('=====================');
save(append(BS2NameResults, '.mat'), ...
    'betaCI', 'sigmaCI', 'tauCI', 'R0BSCI');

cd(currDir)

end