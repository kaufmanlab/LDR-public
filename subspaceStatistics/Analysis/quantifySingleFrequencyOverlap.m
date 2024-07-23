function output = quantifySingleFrequencyOverlap
%% OVERVIEW

% This function is simply the analysis version of 
% quantifyProjectedVarianceSingleFrequency, which mean results are
% successively loaded into the function and returned as results.

%% Analyse data.

% Load the TBF results.
load('TBFResults.mat');

% For each monkey and brain region, analyze.
for monkey = 1:size(TBFResults,2)
    output(monkey).M1 = ...
        quantifyProjectedVarianceSingleFrequency( ...
        TBFResults(monkey).M1.loadings,TBFResults(monkey).M1.basisFxns, ...
        TBFResults(monkey).M1.params.param);
    output(monkey).PMd = ...
        quantifyProjectedVarianceSingleFrequency( ...
        TBFResults(monkey).PMd.loadings,TBFResults(monkey).PMd.basisFxns, ...
        TBFResults(monkey).PMd.params.param);
end

end