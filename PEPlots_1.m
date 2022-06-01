function PEPlots_1(tPeriods, timeVect, dataDate, inputData, inputEst, ...
    outputCurves, outputIncid, outputRepNo, minG, areaPop, LHSNameResults, ...
    days2test)

currDir = cd;
figsDir = 'lhs_figs_1';
cd(figsDir)
if ~exist(LHSNameResults, 'dir')
   mkdir(LHSNameResults)
end
cd(LHSNameResults)

g = minG;

color_model = 1/255*[30,144,255];
color_cases = 1/255*[255,160,122];
color_casesra = 1/255*[255,127,80];
color_S = 1/255*[255,218,185];
color_E = 1/255*[218,165,32];
color_Ir = 1/255*[255,160,122];
color_Iu = 1/255*[255,69,0]';
color_R = 1/255*[144,238,144];
color_V = 1/255*[5,19,171];                                 
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


% Daily data vs fit
figure(1)
tiledlayout(1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
bar(timeVect, inputData(:, 1), 'FaceColor', color_cases, ...
    'EdgeColor', 'none')
hold on
plot(timeVect, inputEst(:, 1), 'color', color_casesra, 'LineWidth', 4);
hold on

CdiffWhole = zeros(1, length(timeVect));
for i = 1:length(tPeriods)
    CdiffWhole(tPeriods{i} + 1) = outputIncid{g, i}(:, 1);
    xline(tPeriods{i}(1), 'color', color_death, 'LineStyle', '--', ...
        'LineWidth', 2)
    hold on
end
plot(timeVect, CdiffWhole, 'color', color_model, 'LineWidth', 4);
hold on

ylabel('new confirmed cases', 'FontSize', 18)
xticks([])
set(gca, 'FontSize', 15)

saveas(gcf, 'PEDaily.fig')
exportgraphics(gcf, 'PEDaily.png', 'Resolution', 300)

% Cumulative data vs fit
figure(2)
tiledlayout(1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
bar(timeVect, cumsum(inputData(:, 1)), 'FaceColor', color_cases, ...
    'EdgeColor', 'none')
hold on
plot(timeVect, cumsum(inputEst(:, 1)), 'color', color_casesra, 'LineWidth', 4);
hold on

CdiffWhole = zeros(1, length(timeVect));
for i = 1:length(tPeriods)
    CdiffWhole(tPeriods{i} + 1) = outputIncid{g, i}(:, 1);
    xline(tPeriods{i}(1), 'color', color_death, 'LineStyle', '--', ...
        'LineWidth', 2)
    hold on
end
plot(timeVect, cumsum(CdiffWhole), 'color', color_model, 'LineWidth', 4);
hold on

ylabel('cumulative confirmed cases', 'FontSize', 18)
xticks([])
set(gca, 'FontSize', 15)

saveas(gcf, 'PECumulative.fig')
exportgraphics(gcf, 'PECumulative.png', 'Resolution', 300)

% Compartments
figure(3)
tiledlayout(6, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
cmpts = {'S', 'E', 'Ir', 'Iu', 'R', 'V'};
compartments = {'susceptible', 'exposed', 'reported', ...
            'unreported', 'recovered', 'vaccinated'};

for j = 1:length(cmpts)
    nexttile;
    XdiffWhole = zeros(1, length(timeVect));
    for i = 1:length(tPeriods)
        Xsol = outputCurves{g, i}(:, j);    
        XdiffWhole(tPeriods{i} + 1) = Xsol;
        xline(tPeriods{i}(1), 'color', color_death, ...
            'LineStyle', '--', 'LineWidth', 2)
        hold on
    end
    plot(timeVect, XdiffWhole, ...
        'color', eval(append('color_', cmpts{j})), ...
        'LineWidth', 4);
    hold on
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

saveas(gcf, 'PECompartments.fig')
exportgraphics(gcf, 'PECompartments.png', 'Resolution', 300)


% Time-varying reproduction number
figure(4)
tiledlayout(1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
RtdiffWhole = zeros(1, length(timeVect));
for i = 1:length(tPeriods)
    Ssol = outputCurves{g, i}(:, 1);
    Nsol = areaPop - outputCurves{g, i}(:, 8);
    RtdiffWhole(tPeriods{i} + 1) = outputRepNo(g, i) .* (Ssol ./ Nsol);
    xline(tPeriods{i}(1), 'color', color_death, ...
        'LineStyle', '--', 'LineWidth', 2)
    hold on
end
plot(timeVect, RtdiffWhole, ...
    'color', color_boot, ...
    'LineWidth', 4);
hold on
ylabel('time-varying reproduction number', 'FontSize', 18)
xlim([0 tPeriods{i}(end)])

yline(1, 'color', color_death, ...
        'LineStyle', '-', 'LineWidth', 4);

xlabel('dates', 'FontSize', 18)
set(gca, 'XTick', dateX - 1, ...
         'XTickLabel', xDate, ...
         'XTickLabelRotation', 45)
set(gca, 'FontSize', 15)

saveas(gcf, 'PERepNo.fig')
exportgraphics(gcf, 'PERepNo.png', 'Resolution', 300)

cd(currDir)

end