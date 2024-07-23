function output = kinFromTBFNoPotent
%% OVERVIEW

% This function uses neural data to predict kinematics from trial-averaged
% firing rates. This is done both using linear regression and nonlinear
% kernel regression. Prediction is bottle-necked through the loading
% matrices for temporal basis factorization, meaning that the orientation
% and initial state of dynamics is related to the movement performed. This
% is done after removing the potent space.

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
load('inferKinematicsFromLoadingsTrialAveraged');

% For each monkey and brain region, prepare for analysis.
for monkey = 1:size(ShenoyMonkeyData,2)
    % Remove M1's potent space, decode.
    ShenoyMonkeyData(monkey).M1 = ...
        removePotent(ShenoyMonkeyData(monkey).M1, ...
        ShenoyMonkeyData(monkey).EMG, ...
        ShenoyMonkeyData(monkey).Kinematics);
    output(monkey).M1 = prepForAnalysis(ShenoyMonkeyData(monkey).M1, ...
        ShenoyMonkeyData(monkey).Kinematics, ...
        kernelType,kernelScale,nonlinearL2,linearL2,holdOutNum, ...
        basisFxnNum(monkey).M1.maxDim,decodeDim, ...
        inferKinematicsFromLoadingsTrialAveraged(monkey).M1);
    % Remove M1's potent space, decode.
    ShenoyMonkeyData(monkey).PMd = ...
        removePotent(ShenoyMonkeyData(monkey).PMd, ...
        ShenoyMonkeyData(monkey).EMG, ...
        ShenoyMonkeyData(monkey).Kinematics);
    output(monkey).PMd = prepForAnalysis(ShenoyMonkeyData(monkey).PMd, ...
        ShenoyMonkeyData(monkey).Kinematics, ...
        kernelType,kernelScale,nonlinearL2,linearL2,holdOutNum, ...
        basisFxnNum(monkey).PMd.maxDim,decodeDim, ...
        inferKinematicsFromLoadingsTrialAveraged(monkey).PMd);
end

end

function output = prepForAnalysis(data,kin,kernelType,kernelScale,nonLinearL2, ...
    linearL2,holdOutNum,factorNum,decodeDim,compareVar)

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
output.linearPrediction.compareVarExplained = compareVar.linearPrediction.varExplained;
% Assign p-vals, means.
output.linearPrediction.velPVal = ...
    signrank(compareVar.linearPrediction.varExplained.vel.array,output.linearPrediction.varExplained.vel.array);
output.linearPrediction.posPVal = ...
    signrank(compareVar.linearPrediction.varExplained.pos.array,output.linearPrediction.varExplained.pos.array);
output.linearPrediction.noPotentWithPotentPos = ...
    [mean(compareVar.linearPrediction.varExplained.pos.array) ...
    mean(output.linearPrediction.varExplained.pos.array)];
output.linearPrediction.noPotentWithPotentVel = ...
    [mean(compareVar.linearPrediction.varExplained.vel.array) ...
    mean(output.linearPrediction.varExplained.vel.array)];

% Perform nonlinear mapping.
[~,output.nonlinearPrediction.predicted] = ...
    predictKinUsingTBFandDataCondStructureNonlinear(kin,holdOutNum,data, ...
    kernelType,[kernelScale nonLinearL2 decodeDim],factorNum);
% Get the variance explained.
output.nonlinearPrediction.varExplained = getVarExplainedKinematics( ...
    output.data,output.nonlinearPrediction.predicted);
output.nonlinearPrediction.compareVarExplained = compareVar.nonlinearPrediction.varExplained;
% Assign p-vals, means.
output.nonlinearPrediction.velPVal = ...
    signrank(compareVar.nonlinearPrediction.varExplained.vel.array,output.nonlinearPrediction.varExplained.vel.array);
output.nonlinearPrediction.posPVal = ...
    signrank(compareVar.nonlinearPrediction.varExplained.pos.array,output.nonlinearPrediction.varExplained.pos.array);
output.nonlinearPrediction.noPotentWithPotentPos = ...
    [mean(compareVar.nonlinearPrediction.varExplained.pos.array) ...
    mean(output.nonlinearPrediction.varExplained.pos.array)];
output.nonlinearPrediction.noPotentWithPotentVel = ...
    [mean(compareVar.nonlinearPrediction.varExplained.vel.array) ...
    mean(output.nonlinearPrediction.varExplained.vel.array)];

end

