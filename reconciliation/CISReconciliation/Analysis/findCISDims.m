function output = findCISDims
%% OVERVIEW

% This function uses Fisher Discriminant Analysis between fixed points and
% initial states to find the CIS.

%% Find CIS.

% Load data.
load('ShenoyMonkeyData');
ShenoyMonkeyData = ShenoyMonkeyData(1:2);
load('TBFResults');

% Loop over datasets.
for monkey = 1:size(ShenoyMonkeyData,2)
    [output(monkey).M1.projection,output(monkey).M1.dimension, ...
        output(monkey).M1.AUC,output(monkey).M1.pVal, ...
        output(monkey).M1.varExplained] = ...
        identifyCIS(pruneRepeats(ShenoyMonkeyData(monkey).M1),TBFResults(monkey).M1.loadings, ...
        TBFResults(monkey).M1.basisFxns,TBFResults(monkey).M1.params.param,5);
    [output(monkey).PMd.projection,output(monkey).PMd.dimension, ...
        output(monkey).PMd.AUC,output(monkey).PMd.pVal, ...
        output(monkey).PMd.varExplained] = ...
        identifyCIS(pruneRepeats(ShenoyMonkeyData(monkey).PMd),TBFResults(monkey).PMd.loadings, ...
        TBFResults(monkey).PMd.basisFxns,TBFResults(monkey).PMd.params.param,5);
end

end