function LHSmatrix = Model_LHS_1(xminlist,xmaxlist,...
                        nsample,distrib,threshold)

%% The results should be compared to the PRCC results section in
%% Supplementary Material D and Table D.1 for different N (specified by
%% "runs" in the script below

%% Sample size N
runs = nsample;

dummymean = 0;
dummysd = 0;
LHSmatrix = zeros(runs, length(xminlist));
if distrib == 'unif'
    for i = 1:length(xminlist)
        xminparam = xminlist(i);
        xmaxparam = xmaxlist(i);
        LHSmatrix(:, i) = LHS_Call_1(xminparam, dummymean, xmaxparam, dummysd, ...
                                    runs, distrib, threshold);
end

%% Save the workspace
save('Model_LHS_1.mat', 'LHSmatrix');
end