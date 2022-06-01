function FCPlots_1(tPeriods, timeVect, dataDate, inputData, inputEst, ...
    minG, areaPop, BSsamples, days2FC, outputFCIncid, outputFCCurves, ...
    outputFCRepNo, beta4FC, days2test, dataCaseFC, ...
    dataCaseMAFC, BSNameResults, ProvinceName, outputSummary)

currDir = cd;
figsDir = 'fc_figs_1';
cd(figsDir)
if ~exist(BSNameResults, 'dir')
   mkdir(BSNameResults)
end
cd(BSNameResults)

csFCC = cumsum(dataCaseFC);
g = minG;
timeVectFC = timeVect(end) + 1: timeVect(end) + days2FC;
xlimVals = [timeVect(end) - days2FC, timeVect(end) + days2FC];
xlimValsFC = [0, timeVect(end) + days2FC];
ylimValsC = [0, max(dataCaseFC(end-30:end))*10 + 1];
ylimValsCA = [0, max(csFCC(end-30:end))*10 + 1];

FCDate = dataDate(1):(dataDate(end) + days2FC);
FC31 = (dataDate(end) + 1):(dataDate(end) + days2FC);
FCtest = (dataDate(end) + 1):(dataDate(end) + days2test);

color_model = 1/255*[30,144,255];
color_cases = 1/255*[255,160,122];
color_casesra = 1/255*[255,127,80];
color_recov = 1/255*[144,238,144];
color_rmod = 1/255*[0,100,0];
color_S = 1/255*[255,218,185];
color_E = 1/255*[218,165,32];
color_C = 1/255*[255,160,122];
color_Iu = 1/255*[255,238,130];
color_R = 1/255*[144,238,144];
color_V = 1/255*[5,19,171];
color_death = 1/255*[112,128,144];
color_dmod = 1/255*[47,79,79];
color_boot = 1/255*[75,0,130];

startMo = dataDate(end - days2FC);
endMo = dataDate(end) + days(days2FC);

dateXLabels = linspace(startMo, endMo, 6);
dateX = datefind(dateXLabels, FCDate);
xDate = datestr(dateXLabels, 'mmm dd');

startMoFC = dateshift(dataDate(1), 'start', 'month');
endMoFC = dateshift(dataDate(end-days2test), 'start', 'month');
if startMoFC < dataDate(1)
    startMoFC = startMoFC + calmonths(1);
    duraMoFC = split(between(startMoFC, endMoFC, 'months'), 'months');
    dateXLabelsFC = sort([startMoFC + calmonths(0:1:duraMoFC)]);
    dateXFC = datefind(dateXLabelsFC, dataDate);
    xDateFC = datestr(dateXLabelsFC, 'mmm yyyy');
else
    duraMoFC = split(between(startMoFC, endMoFC, 'months'), 'months');
    dateXLabelsFC = sort([startMoFC + calmonths(0:1:duraMoFC)]);
    dateXFC = datefind(dateXLabelsFC, dataDate);
    xDateFC = datestr(dateXLabelsFC, 'mmm yyyy');
end


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
    for i = 1:(length(tPeriods))
        CAll(g, tPeriods{i} + 1) = outputFCIncid{g, i}(:, 1);

        SAll(g, tPeriods{i} + 1) = outputFCCurves{g, i}(:, 1);
        EAll(g, tPeriods{i} + 1) = outputFCCurves{g, i}(:, 2);
        IrAll(g, tPeriods{i} + 1) = outputFCCurves{g, i}(:, 3);
        IuAll(g, tPeriods{i} + 1) = outputFCCurves{g, i}(:, 4);
        RAll(g, tPeriods{i} + 1) = outputFCCurves{g, i}(:, 5);
        VAll(g, tPeriods{i} + 1) = outputFCCurves{g, i}(:, 6);
        DAll(g, tPeriods{i} + 1) = outputFCCurves{g, i}(:, 8);

        Ssol = outputFCCurves{g, i}(:, 1);
        Nsol = areaPop - outputFCCurves{g, i}(:, 8);
        R0All(g, tPeriods{i} + 1) = outputFCRepNo(g, i) .* (Ssol ./ Nsol);
    end
end

j = length(tPeriods);

CAllFC = zeros(BSsamples, days2FC, 3);
SAllFC = zeros(BSsamples, days2FC, 3);
EAllFC = zeros(BSsamples, days2FC, 3);
IrAllFC = zeros(BSsamples, days2FC, 3);
IuAllFC = zeros(BSsamples, days2FC, 3);
RAllFC = zeros(BSsamples, days2FC, 3);
VAllFC = zeros(BSsamples, days2FC, 3);
DAllFC = zeros(BSsamples, days2FC, 3);
R0AllFC = zeros(BSsamples, days2FC, 3);
for g = 1:BSsamples
    for i = 1:3
        CAllFC(g, :, i) = outputFCIncid{g, j + i}(:, 1);

        SAllFC(g, :, i) = outputFCCurves{g, j + i}(:, 1);
        EAllFC(g, :, i) = outputFCCurves{g, j + i}(:, 2);
        IrAllFC(g, :, i) = outputFCCurves{g, j + i}(:, 3);
        IuAllFC(g, :, i) = outputFCCurves{g, j + i}(:, 4);
        RAllFC(g, :, i) = outputFCCurves{g, j + i}(:, 5);
        VAllFC(g, :, i) = outputFCCurves{g, j + i}(:, 6);
        DAllFC(g, :, i) = outputFCCurves{g, j + i}(:, 8);

        SsolFC = outputFCCurves{g, j + i}(:, 1);
        NsolFC = areaPop - outputFCCurves{g, j + i}(:, 8);
        R0AllFC(g, :, i) = outputFCRepNo(g, j + i) .* (SsolFC ./ NsolFC);
    end
end

NTV = [timeVect timeVect(end) + (1:days2test)];

% Daily data vs fit
figure(10)
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

for i = 1:3
    nexttile;
    bar(timeVect, inputData(:, 1), 'FaceColor', color_cases, ...
        'EdgeColor', 'none')
    hold on
    semilogy(timeVect, inputEst(:, 1), 'color', color_casesra, 'LineWidth', 4);
    grid on
    hold on

    bar(timeVect(end) + (1:days2test), dataCaseFC(end-days2test+1:end), 'FaceColor', 'none', ...
        'EdgeColor', color_cases)
    hold on
     semilogy(timeVect(end) + (1:days2test), dataCaseMAFC(end-days2test+1:end), ...
         'color', color_casesra, 'LineWidth', 4, 'LineStyle', '--');
    grid on
    hold on
    
    cd(currDir)
    CMed = plims_1(CAll, 0.5);
    CUp = plims_1(CAll, 0.975);
    CLow = plims_1(CAll, 0.025);
    cd(figsDir)
    cd(BSNameResults)   

    casesFit = semilogy(timeVect, CMed, '-', 'Color', color_model, ...
        'LineWidth', 2);
    hold on
    
    fillX = [timeVect, fliplr(timeVect)];
    fillY = [CLow, fliplr(CUp)];
    h = fill(fillX, fillY, color_model, 'linestyle', 'none');
    set(h,'facealpha', .5)
    hold on
    
    cd(currDir)
    CMedFC = plims_1(CAllFC(:, :, i), 0.5);
    CUpFC = plims_1(CAllFC(:, :, i), 0.975);
    CLowFC = plims_1(CAllFC(:, :, i), 0.025);
    cd(figsDir)
    cd(BSNameResults)

    if i == 1
        C4UIT(:, i) = [CMedFC];
        C4UIT(:, i + 1) = [CLowFC];
        C4UIT(:, i + 2) = [CUpFC];
    else
        C4UIT(:, 3*i - 2) = [CMedFC];
        C4UIT(:, 3*i - 1) = [CLowFC];
        C4UIT(:, 3*i) = [CUpFC];
    end

    casesFC = semilogy(timeVectFC, CMedFC, '-', 'Color', color_model, ...
    'LineWidth', 2, 'LineStyle', '--');
    hold on
    
    fillX = [timeVectFC, fliplr(timeVectFC)];
    fillY = [CLowFC, fliplr(CUpFC)];
    h = fill(fillX, fillY, color_model, 'linestyle', 'none');
    set(h,'facealpha', .25)
    hold on
    
    for l = 1:length(tPeriods)
        xline(tPeriods{l}(1), 'color', color_death, 'LineStyle', '--', ...
            'LineWidth', 2)
        hold on
    end
    xline(timeVect(end) + 1, 'color', color_death, 'LineStyle', '-', ...
    'LineWidth', 2)
    
    if i == 1
        ylabel('new confirmed cases', 'FontSize', 18)
    end

    xlim(xlimVals)
    xticks([])
    set(gca, 'FontSize', 15)

    title(append('\beta = ', num2str(beta4FC(i))));
end

saveas(gcf, 'FCDaily.fig')
exportgraphics(gcf, 'FCDaily.png', 'Resolution', 300)


% Cumulative data vs fit
figure(11)
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

for i = 1:3
    nexttile;
    cs = cumsum(inputData(:, 1));
    bar(timeVect, cs, 'FaceColor', color_cases, ...
        'EdgeColor', 'none')
    hold on
    plot(timeVect, cumsum(inputEst(:, 1)), 'color', color_casesra, 'LineWidth', 4);
    grid on
    hold on
    
    csFC = cumsum(dataCaseFC);
    csMAFC = cumsum(dataCaseMAFC);
    bar(timeVect(end) + (1:days2test), csFC(end-days2test+1:end), 'FaceColor', 'none', ...
        'EdgeColor', color_cases)
    hold on
    semilogy(timeVect(end) + (1:days2test), csMAFC(end-days2test+1:end), ...
        'color', color_casesra, 'LineWidth', 4, 'LineStyle', '--');
    grid on
    hold on

    csC = cumsum(CAll, 2);
    cd(currDir)
    CMed = plims_1(csC, 0.5);
    CUp = plims_1(csC, 0.975);
    CLow = plims_1(csC, 0.025);
    cd(figsDir)
    cd(BSNameResults)

    casesFit = semilogy(timeVect, CMed, '-', 'Color', color_model, ...
        'LineWidth', 2);
    hold on
    
    fillX = [timeVect, fliplr(timeVect)];
    fillY = [CLow, fliplr(CUp)];
    h = fill(fillX, fillY, color_model, 'linestyle', 'none');
    set(h,'facealpha', .5)
    hold on
    
    cd(currDir)
    CMedFC = plims_1(csC(:, end) + cumsum(CAllFC(:, :, i), 2), 0.5);
    CUpFC = plims_1(csC(:, end) + cumsum(CAllFC(:, :, i), 2), 0.975);
    CLowFC = plims_1(csC(:, end) + cumsum(CAllFC(:, :, i), 2), 0.025);
    cd(figsDir)
    cd(BSNameResults)

    casesFC = semilogy(timeVectFC, CMedFC, '-', 'Color', color_model, ...
    'LineWidth', 2, 'LineStyle', '--');
    hold on
    
    fillX = [timeVectFC, fliplr(timeVectFC)];
    fillY = [CLowFC, fliplr(CUpFC)];
    h = fill(fillX, fillY, color_model, 'linestyle', 'none');
    set(h,'facealpha', .25)
    hold on
    
    for l = 1:length(tPeriods)
        xline(tPeriods{l}(1), 'color', color_death, 'LineStyle', '--', ...
            'LineWidth', 2)
        hold on
    end
    xline(timeVect(end) + 1, 'color', color_death, 'LineStyle', '-', ...
    'LineWidth', 2)
    
    if i == 1
        ylabel('new confirmed cases', 'FontSize', 18)
    end
    xlim(xlimVals)
    xticks([])
    set(gca, 'FontSize', 15)

    title(append('\beta = ', num2str(beta4FC(i))));
end

saveas(gcf, 'FCCumulative.fig')
exportgraphics(gcf, 'FCCumulative.png', 'Resolution', 300)


% Compartments
figure(12)
tiledlayout(6, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
cmpts = {'S', 'E', 'C', 'Iu', 'R', 'V'};
cmpts1 = {'S', 'E', 'Ir', 'Iu', 'R', 'V'};
compartments = {'susceptible', 'exposed', 'reported', ...
            'unreported', 'recovered', 'vaccinated'};

for j = 1:length(cmpts)
    for i = 1:3
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
        h = fill(fillX, fillY, eval(append('color_', cmpts{j})), 'linestyle', 'none');
        set(h,'facealpha', .5)
        hold on
        
        XVals = eval(append(cmpts1{j}, 'AllFC'));

        cd(currDir)
        XMedFC = plims_1(XVals(:, :, i), 0.5);
        XUpFC = plims_1(XVals(:, :, i), 0.975);
        XLowFC = plims_1(XVals(:, :, i), 0.025);
        cd(figsDir)
        cd(BSNameResults)

        casesFC = plot(timeVectFC, XMedFC, '-', 'Color', eval(append('color_', cmpts{j})), ...
        'LineWidth', 2, 'LineStyle', '--');
        hold on
        
        fillX = [timeVectFC, fliplr(timeVectFC)];
        fillY = [XLowFC, fliplr(XUpFC)];
        h = fill(fillX, fillY, eval(append('color_', cmpts{j})), 'linestyle', 'none');
        set(h,'facealpha', .25)
        hold on

        for l = 1:length(tPeriods)
            xline(tPeriods{l}(1), 'color', color_death, 'LineStyle', '--', ...
                'LineWidth', 2)
            hold on
        end
        xline(timeVect(end) + 1, 'color', color_death, 'LineStyle', '-', ...
        'LineWidth', 2)
    
        if j == 1
            title(append('\beta = ', num2str(beta4FC(i))));
        end

        if i == 1
            ylabel(compartments{j}, 'FontSize', 18)
        end
        xlim(xlimVals)
        set(gca, 'FontSize', 15)
        set(gca, 'XTick', [])

        if i == 2 && j == length(cmpts)
            xlabel('dates', 'FontSize', 18)
        end

        if j == length(cmpts)
            set(gca, 'XTick', dateX - 1, ...
                     'XTickLabel', xDate, ...
                     'XTickLabelRotation', 45)
        else
            xticks([])
        end

        set(gca, 'FontSize', 15)
    end
end

saveas(gcf, 'FCCompartments.fig')
exportgraphics(gcf, 'FCCompartments.png', 'Resolution', 300)



% Time-varying reproduction number
figure(13)
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

for i = 1:3
    nexttile;

    cd(currDir)
    R0Med = plims_1(R0All, 0.5);
    R0Up = plims_1(R0All, 0.975);
    R0Low = plims_1(R0All, 0.025);
    cd(figsDir)
    cd(BSNameResults)

    casesFit = semilogy(timeVect, R0Med, '-', 'Color', color_model, ...
        'LineWidth', 2);
    hold on
    
    fillX = [timeVect, fliplr(timeVect)];
    fillY = [R0Low, fliplr(R0Up)];
    h = fill(fillX, fillY, color_boot, 'linestyle', 'none');
    set(h,'facealpha', .5)
    hold on
    
    cd(currDir)
    R0MedFC = plims_1(R0AllFC(:, :, i), 0.5);
    R0UpFC = plims_1(R0AllFC(:, :, i), 0.975);
    R0LowFC = plims_1(R0AllFC(:, :, i), 0.025);
    cd(figsDir)
    cd(BSNameResults)

    casesFC = semilogy(timeVectFC, R0MedFC, '-', 'Color', color_boot, ...
    'LineWidth', 2, 'LineStyle', '--');
    grid on
    hold on
    
    fillX = [timeVectFC, fliplr(timeVectFC)];
    fillY = [R0LowFC, fliplr(R0UpFC)];
    h = fill(fillX, fillY, color_boot, 'linestyle', 'none');
    set(h,'facealpha', .25)
    hold on

    for l = 1:length(tPeriods)
        xline(tPeriods{l}(1), 'color', color_death, 'LineStyle', '--', ...
            'LineWidth', 2)
        hold on
    end
    xline(timeVect(end) + 1, 'color', color_death, 'LineStyle', '-', ...
    'LineWidth', 2)

    if i == 1
        ylabel('time-varying reproduction number', 'FontSize', 18)
    end
    
    yline(1, 'color', color_death, ...
            'LineStyle', '-', 'LineWidth', 4);
    
    if i == 2
        xlabel('dates', 'FontSize', 18)
    end
    set(gca, 'XTick', dateX - 1, ...
             'XTickLabel', xDate, ...
             'XTickLabelRotation', 45)
    set(gca, 'FontSize', 15)

        R0UIT(:, 3*i - 2) = [R0MedFC];
        R0UIT(:, 3*i - 1) = [R0LowFC];
        R0UIT(:, 3*i) = [R0UpFC];
    title(append('\beta = ', num2str(beta4FC(i))));
    xlim(xlimVals)
%     ylim(ylimValsRo)
end

saveas(gcf, 'FCRepNo.fig')
exportgraphics(gcf, 'FCRepNo.png', 'Resolution', 300)



figure(14)
tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

for i = 1:3
    nexttile;
    bar(timeVect, inputData(:, 1), 'FaceColor', color_cases, ...
        'EdgeColor', 'none')
    hold on
    semilogy(timeVect, inputEst(:, 1), 'color', color_casesra, 'LineWidth', 4);
    grid on
    hold on

    bar(timeVect(end) + (1:days2test), dataCaseFC(end-days2test+1:end), 'FaceColor', 'none', ...
        'EdgeColor', color_cases)
    hold on
    semilogy(timeVect(end) + (1:days2test), dataCaseMAFC(end-days2test+1:end), ...
         'color', color_casesra, 'LineWidth', 4, 'LineStyle', '--');
    grid on
    hold on
    
    cd(currDir)
    CMed = plims_1(CAll, 0.5);
    CUp = plims_1(CAll, 0.975);
    CLow = plims_1(CAll, 0.025);
    cd(figsDir)
    cd(BSNameResults)   

    casesFit = semilogy(timeVect, CMed, '-', 'Color', color_model, ...
        'LineWidth', 2);
    hold on
    
    fillX = [timeVect, fliplr(timeVect)];
    fillY = [CLow, fliplr(CUp)];
    h = fill(fillX, fillY, color_model, 'linestyle', 'none');
    set(h,'facealpha', .5)
    hold on
    
    cd(currDir)
    CMedFC = plims_1(CAllFC(:, :, i), 0.5);
    CUpFC = plims_1(CAllFC(:, :, i), 0.975);
    CLowFC = plims_1(CAllFC(:, :, i), 0.025);
    cd(figsDir)
    cd(BSNameResults)

    casesFC = semilogy(timeVectFC, CMedFC, '-', 'Color', color_model, ...
    'LineWidth', 2, 'LineStyle', '--');
    hold on
    
    fillX = [timeVectFC, fliplr(timeVectFC)];
    fillY = [CLowFC, fliplr(CUpFC)];
    h = fill(fillX, fillY, color_model, 'linestyle', 'none');
    set(h,'facealpha', .25)
    hold on
    
    for l = 1:length(tPeriods)
        xline(tPeriods{l}(1), 'color', color_death, 'LineStyle', '--', ...
            'LineWidth', 2)
        hold on
    end
    xline(timeVect(end) + 1, 'color', color_death, 'LineStyle', '-', ...
    'LineWidth', 2)
    
    if i == 1
        ylabel('new confirmed cases', 'FontSize', 18)
    end
    ylim(ylimValsC)
    xlim(xlimValsFC)
    xticks([])
    set(gca, 'FontSize', 15)

    title(append('\beta = ', num2str(beta4FC(i))));


    ylabel('new confirmed cases', 'FontSize', 18)
    if i == 3
        xlabel('dates', 'FontSize', 18)
    

    set(gca, 'XTick', dateXFC - 1, ...
             'XTickLabel', xDateFC, ...
             'XTickLabelRotation', 45)
    set(gca, 'FontSize', 15)
    xlim(xlimValsFC)
    end
end
saveas(gcf, 'FCDailyextended.fig')
exportgraphics(gcf, 'FCDailyextended.png', 'Resolution', 300)

cd(currDir);

ActualCases = transpose(dataCaseFC(end-days2test+1:end));

FCNameResults = append('Result_FC_',ProvinceName,'_',string(BSsamples),'_',...
                            string(today('datetime')), '_', ...
                             outputSummary);

disp('File name of results:');
disp(FCNameResults);
disp('=====================');

save(append(FCNameResults, '.mat'), ...
    'C4UIT', 'R0UIT', 'ActualCases');

cd(currDir)

end