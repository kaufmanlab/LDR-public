function output = quantifyTBFStats
%% OVERVIEW

% This function quantifies the performance of TBF on neural data. In
% particular, this function quantifies the variance explained in the data,
% the recoverability of the temporal basis functions from the data, and
% the correlations of the previous two metrics with reach duration. 

% There are then compared to two classes of null distributions, the first
% where trials have been scrambled in time to preserve spike-correlations
% and total spike counts but not temporal characteristics, and one where
% trials have been warped to all be the same length. 

%% Quantify statistics.

% Load the data.
load('ShenoyMonkeyData');
ShenoyMonkeyData = ShenoyMonkeyData(1:2);
load('ShenoyMonkeyDataSingleTrial');

% Load the warped datasets.
load('warpedData');
load('warpedDataOtherMethod');

% Load the homogenous data.
load('scrambledData');

% Load the number of basis functions.
load('basisFxnNum');

% For each monkey, analyze.
for monkey = 1:size(ShenoyMonkeyData,2)
    output(monkey).M1 = assignAnalysis( ...
        pruneRepeats(ShenoyMonkeyData(monkey).M1), ...
        warpedData(monkey).M1,warpedDataOtherMethod(monkey).M1, ...
        scrambledData(monkey).M1, ...
        basisFxnNum(monkey).M1.maxDim);
    output(monkey).PMd = assignAnalysis( ...
        pruneRepeats(ShenoyMonkeyData(monkey).PMd), ...
        warpedData(monkey).PMd,warpedDataOtherMethod(monkey).PMd, ...
        scrambledData(monkey).PMd, ...
        basisFxnNum(monkey).PMd.maxDim);
end

end

%% FUNCTION FOR ASSIGNING ANALYSIS.

function output = assignAnalysis(data,warped,warped2,scrambled,dim)
% Assign the durations.
output.durations = [warped().oldDuration];
% Analyze the data.
output.data = analysis(data,dim);
output.warped = analysis(warped,dim);
output.warped2 = analysis(warped2,dim);
output.scrambled = analysis(scrambled,dim);
% Get the pvals for scrambled. vs data for variance explained.
output.nullMeanVarExplained = mean(output.scrambled.varExplained);
output.nullStdVarExplained = std(output.scrambled.varExplained);
output.dataMeanVarExplained = mean(output.data.varExplained);
output.dataStdVarExplained = std(output.data.varExplained);
output.varExplainedPVal = ranksum( ...
    output.data.varExplained,output.scrambled.varExplained);
% Get the pvals for scrambled. vs data for recoverability.
output.nullMeanRecoverability = mean(output.scrambled.recoverability);
output.nullStdRecoverability = std(output.scrambled.recoverability);
output.dataMeanRecoverability = mean(output.data.recoverability);
output.dataStdRecoverability = std(output.data.recoverability);
output.recoverabilityPVal = ranksum( ...
    output.data.recoverability,output.scrambled.recoverability);
% Get the correlation between recoverability and the duration in the data
% and the warped with the first method and the second method.
[output.dataRecoverabilityCorr,output.dataRecoverabilityP] = ...
    corr(abs(output.durations.'-mean(output.durations)),output.data.recoverability.');
[output.nullRecoverabilityCorr,output.nullRecoverabilityP] = ...
    corr(abs(output.durations.'-mean(output.durations)),output.warped.recoverability.');
[output.null2RecoverabilityCorr,output.null2RecoverabilityP] = ...
    corr(abs(output.durations.'-mean(output.durations)),output.warped2.recoverability.');
end

%% Function for analyzing data.

function output = analysis(data,dim)
% Get the variance explained and recoverability for the normal data.
[basisFxns,loadings] = getBasisFxns(data,dim);
output.varExplained = getVarExplained(data,useFactorization(loadings,basisFxns),'ind');
output.varExplained = output.varExplained.array;
output.recoverability = getRecoverability( ...
        data,loadings,basisFxns);
output.recoverability = output.recoverability.array;

end










