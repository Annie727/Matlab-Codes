function prccPlot = plotVariedPRCC_1(M,N,x,labelstring,parameters,prcc)

linS = {'-','-','--',':','-','--','-','-',':',':',':','--','-'};
% linS = {'-','-','--',':','-','--'};
    clr = {[0.83 0.14 0.14],...
           [0.47 0.25 0.80],...
           [0.47 0.25 0.80],...
           [1.00 1.00 0.00],...
           [0.7  0.7  0.7],...
           [0.25 0.80 0.54],...
           [0.25 0.80 0.54],...
           [0.00 0.00 0.00],...
           [0.3010 0.7450 0.9330],...
           [1.00 0.078 0.576],...
           [1.00 0.54 0.00],...
           [0.00 0.54 0.7],...
           [0.25 1.00 0.14]};

% clr = {[0.83 0.14 0.14],...
%            [0.47 0.25 0.80],...
%            [0.47 0.25 0.80],...
%            [1.00 1.00 0.00],...
%            [0.7  0.7  0.7],...
%            [0.25 0.80 0.54]};
    figure()
    hold on
    box on
    for mm=1:M
        plot(x, prcc(mm,:),'LineWidth',2.0,'linestyle',linS{mm},'Color',clr{mm});
        
            
    end

    prccPlot = gca; 

    xlabel('time (days)');
    ylabel('PRCC value');
    xlim([0, x(end)])
    title(['Correlation Plot of ',labelstring,' from LHS simulations, ',num2str(N),' samples']);
    legend(parameters.name,'Location','EastOutside','Fontsize',15)
    xdate=datestr(datenum('03/01/2021'):30:datenum('11/30/2021'));
    set(gca,'XTick',0:30:length(x)-1,'XTickLabel',xdate,'XTickLabelRotation',45)
end