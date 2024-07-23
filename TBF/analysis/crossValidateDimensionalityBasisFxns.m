function output = crossValidateDimensionalityBasisFxns
%% OVERVIEW

% This is an analysis function, meaning it loads its own data. This
% function is responsible for cross-validating the dimensionality of the
% temporal basis functions. This is done by partitioning the trials into
% two portions, performing temporal basis factorization on varying
% dimensionalities of the basis functions, and then finding the
% dimensionality that minimizes the empirical loss on the training set. 

% Sorry everything is called 1 I just don't want to edit it out. 

%% Parameters

% Number of folds to use for method 1, the train ratio for trials.
method1Folds = 100;
method1Ratio = 1/2;

%% Cross-validate.

% Load the dataset.
load('ShenoyMonkeyData');
ShenoyMonkeyData = ShenoyMonkeyData(1:2);
load('ShenoyMonkeyDataSingleTrial');

% For each monkey, cross validate each region.
for monkey = 1:size(ShenoyMonkeyData,2)
    if ShenoyMonkeyData(monkey).singleTrialExist
        trialM1 = ShenoyMonkeyDataSingleTrial(monkey).M1;
        trialPMd = ShenoyMonkeyDataSingleTrial(monkey).PMd;
    else
        trialM1 = [];
        trialPMd = [];
    end
    output(monkey).M1 = assignMethod(ShenoyMonkeyData(monkey).M1,trialM1, ...
        [method1Folds method1Ratio]);
    output(monkey).PMd = assignMethod(ShenoyMonkeyData(monkey).PMd,trialPMd, ...
        [method1Folds method1Ratio]);
end

end

%% ASSIGN METHODS

function output = assignMethod(condData,trialData,method1Params)
% Prune the data.
condData = pruneRepeats(condData);
% Cross validate.
output = method1(condData,trialData,method1Params(1),method1Params(2));
end

%% METHOD OF CROSSVALIDATING

function output = method1(data,trialData,method1Folds,method1Ratio)
% Process trials.
if isempty(trialData)
    trials = generateTrials(data);
else
    trials = pruneRepeatsTrial(trialData);
end
% Pre-process the trials.
trials = preprocessTrials(trials);
% Assign the dimensionalities to test.
dimTest = 1:20;
output.lossArray = zeros(length(dimTest),method1Folds);
% Fit dynamics to train conditions.
for fold = 1:method1Folds
    % Partition into train and test conditions.
    [trainData,testData] = cvConds(trials,method1Ratio);
    for dim = dimTest
        [basisFxns,loadings] = getBasisFxns(trainData,dim);
        % Assign the variance explained.
        VE = getVarExplained(testData,useFactorization(loadings,basisFxns),'ind');
        output.lossArray(dim,fold) = ...
            VE.mean;
    end
end
% Quantify stats.
output.meanTrace = ...
    mean(output.lossArray,2);
output.stdTrace = ...
    std(output.lossArray,0,2);
[~,output.maxDim] = max(output.meanTrace);
if ~mod(output.maxDim,2)
    output.maxDim = output.maxDim+1;
end
end

