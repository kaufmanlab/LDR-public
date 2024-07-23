function output = initialStatesToActivity
%% OVERVIEW

% This function attempts to predict the evolution of motor cortex activity
% on a given condition from the initial state using linear regression.
 
%% Parameters.

% Refers to the regularization on the regression.
L2 = 0.1;

% Refers to the portion of the singular value spectrum of the initial
% states to keep.
svEnergy = 0.8;

% Refers to the portion of the condition to leave out from
% cross-validation.
leaveOut = 6;

%% Fit dynamics.

% Load data.
load('ShenoyMonkeyData');
ShenoyMonkeyData = ShenoyMonkeyData(1:2);
load('basisFxnNum');

% Loop over monkeys. 
for monkey = 1:size(ShenoyMonkeyData,2)
    pooled = poolDatasets(ShenoyMonkeyData(monkey).M1,ShenoyMonkeyData(monkey).PMd);
    output(monkey).M1 = analyze(ShenoyMonkeyData(monkey).M1,L2,svEnergy, ...
        leaveOut,basisFxnNum(monkey).M1.maxDim,pooled);
    output(monkey).PMd = analyze(ShenoyMonkeyData(monkey).PMd,L2,svEnergy, ...
        leaveOut,basisFxnNum(monkey).PMd.maxDim,pooled);
end

end

%% SUBFUNCTION FOR ANALYZING DATA.

function out = analyze(data,L2,svEnergy,leaveOut,factorDim,pooled)
data = pruneRepeats(data);
pooled = pruneRepeats(pooled);
% Get the initial states.
initialStates = zeros(size(pooled(1).matrix,1),size(pooled,2));
for cond = 1:size(pooled,2)
    initialStates(:,cond) = mean(pooled(cond).matrix(:,1:10),2);
end
% Predict activity from initial states.
[out.data,out.prediction] = predictActivityFromInitialState( ...
    data,leaveOut,initialStates,L2,factorDim,svEnergy);
% Quantify the variance explained.
out.varExplained = getVarExplained(out.data,out.prediction,'ind');
% Predict activity from shuffled initial states.
[nullData,nullPrediction] = predictActivityFromInitialState( ...
    data,leaveOut,initialStates(:,randperm(size(initialStates,2))), ...
    L2,factorDim,svEnergy);
% Quantify the variance explained.
out.nullVarExplained = getVarExplained(nullData,nullPrediction,'ind');
% Quantify the p value.
out.pVal = signrank(out.nullVarExplained.array,out.varExplained.array);
end

