function output = fitDynamics
%% OVERVIEW

% This function attempts to predict the evolution of motor cortex activity
% on a given condition from previously-observed conditions using NN
% regression. This is done by predicting the derivative of the state from
% the state itself. 

%% Fit dynamics.

% Load data.
load('ShenoyMonkeyData');
ShenoyMonkeyData = ShenoyMonkeyData(1:2);
load('TBFResults');

% Loop over monkeys. 
for monkey = 1:size(ShenoyMonkeyData,2)
    output(monkey).M1 = analyze(ShenoyMonkeyData(monkey).M1, ...
        TBFResults(monkey).M1);
    output(monkey).PMd = analyze(ShenoyMonkeyData(monkey).PMd, ...
        TBFResults(monkey).PMd);
end

end

%% SUBFUNCTION FOR ANALYZING DATA.

function out = analyze(data,resultsTBF)
data = pruneRepeats(data);
% Fit dynamics.
[out.prediction,out.data] = ...
    predictDynamicsUsingNNRegression(data,6,resultsTBF.loadings,resultsTBF.basisFxns);
% Quantify the variance explained.
out.varExplained = getVarExplained(out.data,out.prediction,'ind');
% Fit dynamics to shuffled data
shuffledLoadings = resultsTBF.loadings;
for cond = 1:size(resultsTBF.loadings,2)
    for col = 1:size(shuffledLoadings(cond).matrix,2)
        randCond = randi([1 size(data,2)],1,1);
    	shuffledLoadings(cond).matrix(:,col) = resultsTBF.loadings(randCond).matrix(:,col);
    end
end
[nullPrediction,nullData] = ...
    predictDynamicsUsingNNRegression(data,6,shuffledLoadings,resultsTBF.basisFxns);
% Quantify the variance explained.
out.nullVarExplained = getVarExplained(nullData,nullPrediction,'ind');
% Quantify the p value.
out.pVal = signrank(out.nullVarExplained.array,out.varExplained.array);
end

