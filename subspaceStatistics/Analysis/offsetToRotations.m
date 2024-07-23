function output = offsetToRotations
%% OVERVIEW

% This function uses neural data to predict kinematics from trial-averaged
% firing rates. This is done both using linear regression and nonlinear
% kernel regression. Prediction is bottle-necked through the loading
% matrices for temporal basis factorization, meaning that the orientation
% and initial state of dynamics is related to the movement performed. 

%% Parameters.

% Declare the kernel used for non-linear mapping.
kernelType = 'OU';

% Declares the length scale of the kernel.
kernelScale = 5*10^2;

% Declare the regularization on the kernel for nonlinear mapping.
nonlinearL2 = 0.00001;

% Declare the dimensionality for decoding.
decodeEnergy = 0.9;

% Declare the regularization on the linear mapping.
linearL2 = 0.01;

% Declare the train-test ratio.
holdOutNum = 2;

%% Form the mapping.

% Load the data.
load('ShenoyMonkeyData');
ShenoyMonkeyData = ShenoyMonkeyData(1:2);
load('basisFxnNum');
load('TBFResults');

% For each monkey and brain region, prepare for analysis.
for monkey = 1:size(ShenoyMonkeyData,2)
    % M1.
    output(monkey).M1 = prepForAnalysis(ShenoyMonkeyData(monkey).M1, ...
        kernelType,kernelScale,nonlinearL2,linearL2,holdOutNum, ...
        basisFxnNum(monkey).M1.maxDim,decodeEnergy,TBFResults(monkey).M1);
    % PMd.
    output(monkey).PMd = prepForAnalysis(ShenoyMonkeyData(monkey).PMd, ...
        kernelType,kernelScale,nonlinearL2,linearL2,holdOutNum, ...
        basisFxnNum(monkey).PMd.maxDim,decodeEnergy,TBFResults(monkey).PMd);
end

end

function output = prepForAnalysis(data,kernelType,kernelScale,nonLinearL2, ...
    linearL2,holdOutNum,factorNum,svEnergy,resultsTBF)

% Prune out repeats.
data = pruneRepeats(data);

% Perform linear mapping.
[output.data,output.linearPrediction.predicted] = ...
    predictRotationsFromOffset( ...
    data,holdOutNum,linearL2,factorNum,svEnergy,1);
% Get the variance explained.
output.linearPrediction.varExplained = getVarExplained( ...
    output.data,output.linearPrediction.predicted,'ind');
% Scramble the data.
fakeData = data;
for cond = 1:size(data,2)
    fakeData(cond).matrix = ...
        [resultsTBF.loadings(randi([1 72],1,1)).matrix(:,1) ...
        resultsTBF.loadings(randi([1 72],1,1)).matrix(:,2:end)] ...
        *resultsTBF.params.rawBasisFxns.';
end
% Get the variance explained after scrambling.
[~,null] = ...
    predictRotationsFromOffset( ...
    fakeData,holdOutNum,linearL2,factorNum,svEnergy,1);
output.linearPrediction.nullVarExplained = getVarExplained( ...
    output.data,null,'ind');
output.linearPrediction.pval = signrank( ...
    output.linearPrediction.nullVarExplained.array, ...
    output.linearPrediction.varExplained.array);

% Perform nonlinear mapping.
[output.data,output.nonlinearPrediction.predicted] = ...
    predictRotationsFromOffsetNonlinear( ...
    data,holdOutNum,nonLinearL2,factorNum,svEnergy,kernelType,kernelScale);
% Get the variance explained.
output.nonlinearPrediction.varExplained = getVarExplained( ...
    output.data,output.nonlinearPrediction.predicted,'ind');
% Get the variance explained after scrambling.
[~,null] = ...
predictRotationsFromOffsetNonlinear( ...
    fakeData,holdOutNum,nonLinearL2,factorNum,svEnergy,kernelType,kernelScale);
output.nonlinearPrediction.nullVarExplained = getVarExplained( ...
    output.data,null,'ind');
output.nonlinearPrediction.pval = signrank( ...
    output.nonlinearPrediction.nullVarExplained.array, ...
    output.nonlinearPrediction.varExplained.array);

end

