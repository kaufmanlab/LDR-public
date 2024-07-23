function output = crossValidateDimensionalityConditions
%% OVERVIEW

% This is an analysis function, meaning it loads its own data. This
% function is responsible for cross-validating the dimensionality of the
% conditions to use in naive condition-specific dynamics. This is done in two
% ways. 

% Where single trials are available (array data), single trial are
% partitioned into train and test sets and respectively averaged. Subspaces
% are fit to the train set, and then the quality of the approximation is
% tested on the test set. This is done multiple times to form a loss
% curve. This method is simulated for single-unit data by sampling "new trials"
% via Poisson sampling of averaged firing rates to establish a bound of the
% possible SNR with trial counts.

% Sorry everything is called "method 1", there used to be a method 2 but I
% axed it and don't want to rename everything.

%% Parameters

% Number of folds to use for method 1, the train ratio for trials.
method1Folds = 100;
method1Ratio = 3/4;

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
dimTest = 1:2:21;
% Assign an array for the loss.
if size(data(1).matrix,1) < dimTest(end)
    maxDim = find(dimTest == ...
        (size(data(1).matrix,1) - 1 + mod(size(data(1).matrix,1),2)));
else
    maxDim = length(dimTest);
end
output.lossArray = zeros(maxDim,size(data,2),method1Folds);
% Fit dynamics to train conditions.
for fold = 1:method1Folds
    % Partition into train and test conditions.
    [trainData,testData] = cvConds(trials,method1Ratio);
    for dim = 1:maxDim
        [~,~,~,approx,~,~] = idDynamics(trainData,dimTest(dim));
        % Assign the variance explained.
        VE = getVarExplained(testData,approx,'ind');
        output.lossArray(dim,:,fold) = ...
            VE.array;
    end
end
% Quantify stats.
output.meanTrace = ...
    mean(reshape(output.lossArray,maxDim,size(data,2)*method1Folds),2);
output.stdTrace = ...
    std(reshape(output.lossArray,maxDim,size(data,2)*method1Folds),0,2);
output.meanTraceCond = ...
    mean(output.lossArray,3);
output.stdTraceCond = ...
    std(output.lossArray,0,3);
[~,output.maxDim] = max(output.meanTrace);
output.maxDim = dimTest(output.maxDim);
[~,output.maxDimCond] = max(output.meanTraceCond);
output.maxDimCond = dimTest(output.maxDimCond);
end

