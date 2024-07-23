function output = crossValidateDimensionalityRotations
%% OVERVIEW

% This function cross-validates the dimensionality of each rotation
% independently. This is done by getting seperate estimates of each
% rotation from independent splits of the data, then the dimensionality of
% the subspace of each rotation is swept to maximize the variance explained
% in the alternate estimate of the rotation.

%% Parameters

% Number of folds to use for method 1, the train ratio for trials.
method1Folds = 100;
method1Ratio = 1/2;

%% Cross-validate.

% Load the dataset.
load('ShenoyMonkeyData');
ShenoyMonkeyData = ShenoyMonkeyData(1:2);
load('ShenoyMonkeyDataSingleTrial');
load('basisFxnNum');

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
        [method1Folds method1Ratio],basisFxnNum(monkey).M1.maxDim);
    output(monkey).PMd = assignMethod(ShenoyMonkeyData(monkey).PMd,trialPMd, ...
        [method1Folds method1Ratio],basisFxnNum(monkey).PMd.maxDim);
end

end

%% ASSIGN METHODS

function output = assignMethod(condData,trialData,method1Params,dim)
% Prune the data.
condData = pruneRepeats(condData);
% Cross validate.
output = method1(condData,trialData,method1Params(1),method1Params(2),dim);
end

%% METHOD OF CROSSVALIDATING

function output = method1(data,trialData,method1Folds,method1Ratio,maxDim)
% Process trials.
if isempty(trialData)
    trials = generateTrials(data);
else
    trials = pruneRepeatsTrial(trialData);
end
% Pre-process the trials.
trials = preprocessTrials(trials);
% Assign the dimensionalities to test.
dimTest = 1:50;
output.lossArray = zeros(length(dimTest),method1Folds,(maxDim-1)/2);
% Fit dynamics to train conditions.
for fold = 1:method1Folds
    % Partition into train and test conditions.
    [trainData,testData] = cvConds(trials,method1Ratio);
    [basisFxns,loadings,params] = eigTransform(trainData,maxDim);
    % Fit test loadings.
    for cond = 1:size(loadings,2)
        loadings(cond).test = testData(cond).matrix*pinv(params.rawBasisFxns.');
    end
    % Get the frequencies present in the fxns.
    uniqueFreqs = unique(params.param(:,1));
    % Loop over rotations
    for rot = 1:length(uniqueFreqs)
        % Get the indices of the rotation.
        findInds = find(params.param(:,1) == uniqueFreqs(rot));
        % Assign an array for the respective columns. 
        colArray = zeros(size(loadings(1).matrix,1),size(loadings,2)*numel(findInds));
        % Populate the array.
        for cond = 1:size(loadings,2)
            colArray(:,(1:numel(findInds))+(cond-1)*numel(findInds)) = loadings(cond).matrix(:,findInds);
        end
        % Get the top singular vectors of the array, to estimate the
        % subspace of the rotations.
        [subspace,~,~] = svd(colArray,'econ');
        % Test the dimensionality.
        for dim = dimTest
            % Assign an array for the variance explained.
            VEArray = zeros(1,size(loadings,2));
            % Loop over conditions.
            for cond = 1:size(loadings,2)
                VEArray(cond) = 1 - sum(var( ...
                    ((subspace(:,1:dim)*subspace(:,1:dim).') ...
                    *loadings(cond).matrix(:,findInds) - ...
                    loadings(cond).test(:,findInds)) ...
                    *basisFxns(21:end,findInds).',0,2)) ...
                    /sum(var(loadings(cond).test(:,findInds) ...
                    *basisFxns(21:end,findInds).',0,2));
            end
            output.lossArray(dim,fold,rot) = ...
                mean(VEArray);
        end
    end
end
% Quantify stats.
output.meanTrace = ...
    squeeze(mean(output.lossArray,2));
output.stdTrace = ...
    squeeze(std(output.lossArray,0,2));
[~,output.maxDim] = max(output.meanTrace);
end

