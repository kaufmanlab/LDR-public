function output = nonlinearControls
%% OVERVIEW

% This function shows that the success of nonlinear mapping is not due to
% memorization of the dataset, but rather due learning a mapping that is
% locally valid. This is done by comparing the performance of the
% regression against nearest neighbors regression, which is by definition
% simple memorization of the training set.  

%% Parameters.

% Declare the kernel used for non-linear mapping.
kernelType = 'OU';

% Declares the length scale of the kernel for [position, velocity].
kernelScale = [10^12 10^11];

% Declare the regularization on the kernel for nonlinear mapping.
nonlinearL2 = 0;

% Declare the regularization on the linear mapping.
linearL2 = 0.1;

% Declare the train-test ratio.
holdOutNum = 6;

%% Form the mapping.

% Load the data.
load('ShenoyMonkeyData');
ShenoyMonkeyData = ShenoyMonkeyData(1:2);
load('inferLoadingsFromKinematicsTrialAveraged');
load('basisFxnNum');

% For each monkey and brain region, prepare for analysis.
for monkey = 1:size(ShenoyMonkeyData,2)
    % M1.
    output(monkey).M1 = prepForAnalysis(ShenoyMonkeyData(monkey).M1, ...
        ShenoyMonkeyData(monkey).Kinematics,holdOutNum, ...
        inferLoadingsFromKinematicsTrialAveraged(monkey).M1.nonlinearPrediction, ...
        basisFxnNum(monkey).M1.maxDim, ...
        kernelType,kernelScale,nonlinearL2,linearL2);
    % PMd.
    output(monkey).PMd = prepForAnalysis(ShenoyMonkeyData(monkey).PMd, ...
        ShenoyMonkeyData(monkey).Kinematics,holdOutNum, ...
        inferLoadingsFromKinematicsTrialAveraged(monkey).PMd.nonlinearPrediction, ...
        basisFxnNum(monkey).PMd.maxDim, ...
        kernelType,kernelScale,nonlinearL2,linearL2);
end

end

%% Subfunction for organizing datasets before analysis.

function output = prepForAnalysis(data,kin,holdOutNum,nonlinearPred,factorNum, ...
    kernelType,kernelScale,nonLinearL2,linearL2)

% Prune out repeats.
data = pruneRepeats(data);
kin = pruneRepeats(kin);

% Assign the nonlinear prediction.
output.nonlinearPrediction = nonlinearPred;

% Perform NN prediction.
[~,output.NN.predicted] = predictDataUsingTBFandKinNN(data,holdOutNum,kin,factorNum,'angles');
output.NN.varExplained = getVarExplained( ...
    data,output.NN.predicted,'ind');
output.NN.pVal = signrank( ...
    output.NN.varExplained.array, ...
    output.nonlinearPrediction.control.varExplained.array);

% Perform linear encoding.
[lagged,output.linearEncoding.predicted] = predictDataUsingKin(data,holdOutNum,kin,linearL2, ...
    10,'rand');
output.linearEncoding.varExplained = getVarExplained( ...
    lagged,output.linearEncoding.predicted,'ind');
output.linearEncoding.pVal = signrank( ...
    output.linearEncoding.varExplained.array, ...
    output.nonlinearPrediction.varExplained.array);

% Perform nonlinear position encoding.
[lagged,output.nonlinearPosition.predicted] = ...
    predictDataUsingKinNonlinear(data,holdOutNum,kin,'Pos', ...
    kernelType,[kernelScale(1) nonLinearL2],10,'rand');
output.nonlinearPosition.varExplained = getVarExplained( ...
    lagged,output.nonlinearPosition.predicted,'ind');
output.nonlinearPosition.pVal = signrank( ...
    output.nonlinearPosition.varExplained.array, ...
    output.nonlinearPrediction.varExplained.array);

% Perform nonlinear velocity encoding.
[lagged,output.nonlinearVelocity.predicted] = ...
    predictDataUsingKinNonlinear(data,holdOutNum,kin,'Vel', ...
    kernelType,[kernelScale(2) nonLinearL2],10,'rand');
output.nonlinearVelocity.varExplained = getVarExplained( ...
    lagged,output.nonlinearVelocity.predicted,'ind');
output.nonlinearVelocity.pVal = signrank( ...
    output.nonlinearVelocity.varExplained.array, ...
    output.nonlinearPrediction.varExplained.array);

end

