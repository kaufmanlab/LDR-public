function output = linearControls
%% OVERVIEW

% This function shows that the success of linear mapping is not due to some
% confound such as linear encoding of kinematics, in which case
% sequence-to-sequence mapping would work but trivially, as the sequences
% would be some linear transform of each other. In this case the linear
% encoding model is rather broad in the parameters included. This function
% also tests against the possibility that of non-linear encoding of
% position, velocity, or complex future-tuning. 

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

% For each monkey and brain region, prepare for analysis.
for monkey = 1:size(ShenoyMonkeyData,2)
    % M1.
    output(monkey).M1 = prepForAnalysis(ShenoyMonkeyData(monkey).M1, ...
        ShenoyMonkeyData(monkey).Kinematics, ...
        kernelType,kernelScale,nonlinearL2,linearL2,holdOutNum, ...
        inferLoadingsFromKinematicsTrialAveraged(monkey).M1.linearPrediction);
    % PMd.
    output(monkey).PMd = prepForAnalysis(ShenoyMonkeyData(monkey).PMd, ...
        ShenoyMonkeyData(monkey).Kinematics, ...
        kernelType,kernelScale,nonlinearL2,linearL2,holdOutNum, ...
        inferLoadingsFromKinematicsTrialAveraged(monkey).PMd.linearPrediction);
end

end

%% Subfunction for organizing datasets before analysis.

function output = prepForAnalysis(data,kin,kernelType,kernelScale,nonLinearL2, ...
    linearL2,holdOutNum,linearPreds)

% Prune out repeats.
data = pruneRepeats(data);
kin = pruneRepeats(kin);

% Assign the linear prediction.
output.linearPrediction = linearPreds;

% Perform linear encoding.
[lagged,output.linearEncoding.predicted] = predictDataUsingKin(data,holdOutNum,kin,linearL2, ...
    10,'rand');
output.linearEncoding.varExplained = getVarExplained( ...
    lagged,output.linearEncoding.predicted,'ind');
output.linearEncoding.pVal = signrank( ...
    output.linearEncoding.varExplained.array, ...
    output.linearPrediction.varExplained.array);

% Perform nonlinear position encoding.
[lagged,output.nonlinearPosition.predicted] = ...
    predictDataUsingKinNonlinear(data,holdOutNum,kin,'Pos', ...
    kernelType,[kernelScale(1) nonLinearL2],10,'rand');
output.nonlinearPosition.varExplained = getVarExplained( ...
    lagged,output.nonlinearPosition.predicted,'ind');
output.nonlinearPosition.pVal = signrank( ...
    output.nonlinearPosition.varExplained.array, ...
    output.linearPrediction.varExplained.array);

% Perform nonlinear velocity encoding.
[lagged,output.nonlinearVelocity.predicted] = ...
    predictDataUsingKinNonlinear(data,holdOutNum,kin,'Vel', ...
    kernelType,[kernelScale(2) nonLinearL2],10,'rand');
output.nonlinearVelocity.varExplained = getVarExplained( ...
    lagged,output.nonlinearVelocity.predicted,'ind');
output.nonlinearVelocity.pVal = signrank( ...
    output.nonlinearVelocity.varExplained.array, ...
    output.linearPrediction.varExplained.array);

end

