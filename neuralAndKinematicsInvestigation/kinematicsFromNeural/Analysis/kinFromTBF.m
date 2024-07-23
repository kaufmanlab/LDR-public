function output = kinFromTBF
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
kernelScale = 10^3;

% Declare the regularization on the kernel for nonlinear mapping.
nonlinearL2 = 0;

% Declare the dimensionality for decoding.
decodeDim = 20;

% Declare the regularization on the linear mapping.
linearL2 = 1;

% Declare the train-test ratio.
holdOutNum = 2;

%% Form the mapping.

% Load the data.
load('ShenoyMonkeyData');
ShenoyMonkeyData = ShenoyMonkeyData(1:2);
load('basisFxnNum');

% For each monkey and brain region, prepare for analysis.
for monkey = 1:size(ShenoyMonkeyData,2)
    % M1.
    output(monkey).M1 = prepForAnalysis(ShenoyMonkeyData(monkey).M1, ...
        ShenoyMonkeyData(monkey).Kinematics, ...
        kernelType,kernelScale,nonlinearL2,linearL2,holdOutNum, ...
        basisFxnNum(monkey).M1.maxDim,decodeDim);
    % PMd.
    output(monkey).PMd = prepForAnalysis(ShenoyMonkeyData(monkey).PMd, ...
        ShenoyMonkeyData(monkey).Kinematics, ...
        kernelType,kernelScale,nonlinearL2,linearL2,holdOutNum, ...
        basisFxnNum(monkey).PMd.maxDim,decodeDim);
end

end

function output = prepForAnalysis(data,kin,kernelType,kernelScale,nonLinearL2, ...
    linearL2,holdOutNum,factorNum,decodeDim)

% Prune out repeats.
data = pruneRepeats(data);
kin = pruneRepeats(kin);

% Perform linear mapping.
[output.data,output.linearPrediction.predicted] = ...
    predictKinUsingTBFandDataCondStructure(kin,holdOutNum, ...
    data,linearL2,factorNum);
% Get the variance explained.
output.linearPrediction.varExplained = getVarExplainedKinematics( ...
    output.data,output.linearPrediction.predicted);
% Get the variance explained after scrambling.
[~,null] = ...
    predictKinUsingTBFandDataCondStructure(kin,holdOutNum, ...
    data(randperm(size(data,2))),linearL2,factorNum);
output.linearPrediction.nullVarExplained = getVarExplainedKinematics( ...
    output.data,null);
output.linearPrediction.pValPos = signrank( ...
    output.linearPrediction.nullVarExplained.pos.array, ...
    output.linearPrediction.varExplained.pos.array);
output.linearPrediction.pValVel = signrank( ...
    output.linearPrediction.nullVarExplained.vel.array, ...
    output.linearPrediction.varExplained.vel.array);

% Perform nonlinear mapping.
[~,output.nonlinearPrediction.predicted] = ...
    predictKinUsingTBFandDataCondStructureNonlinear(kin,holdOutNum,data, ...
    kernelType,[kernelScale nonLinearL2 decodeDim],factorNum);
% Get the variance explained.
output.nonlinearPrediction.varExplained = getVarExplainedKinematics( ...
    output.data,output.nonlinearPrediction.predicted);
% Get the variance explained after scrambling.
[~,null] = ...
    predictKinUsingTBFandDataCondStructureNonlinear(kin,holdOutNum, ...
    data(randperm(size(data,2))), ...
    kernelType,[kernelScale nonLinearL2 decodeDim],factorNum);
output.nonlinearPrediction.nullVarExplained = getVarExplainedKinematics( ...
    output.data,null);
output.nonlinearPrediction.pValPos = signrank( ...
    output.nonlinearPrediction.nullVarExplained.pos.array, ...
    output.nonlinearPrediction.varExplained.pos.array);
output.nonlinearPrediction.pValVel = signrank( ...
    output.nonlinearPrediction.nullVarExplained.vel.array, ...
    output.nonlinearPrediction.varExplained.vel.array);

end

