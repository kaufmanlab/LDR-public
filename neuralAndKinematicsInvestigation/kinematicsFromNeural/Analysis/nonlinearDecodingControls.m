function output = nonlinearDecodingControls
%% OVERVIEW

% This function controls for the nonlinear decoding performed using TBF by
% comparing it linear decoding from recordings, and non-linear decoding from 
% recordings. Instantaneous decoding is compared with both 20 ms and 100 ms smoothing. 

%% Parameters.

% Declare the kernel used for non-linear mapping.
kernelType = 'OU';

% Declares the length scale of the kernel.
kernelScale = 10^2;

% Declare the regularization on the kernel for nonlinear mapping.
nonlinearL2 = 0.1;

% Declare the regularization on the linear mapping.
linearL2 = 1;

% Declare the train-test ratio.
holdOutNum = 2;

%% Control for decoding.

% Load the data.
load('ShenoyMonkeyData');
ShenoyMonkeyData = ShenoyMonkeyData(1:2);
load('inferKinematicsFromLoadingsTrialAveraged');

% For each monkey and brain region, prepare for analysis.
for monkey = 1:size(ShenoyMonkeyData,2)
    % M1.
    output(monkey).M1 = prepForAnalysis(ShenoyMonkeyData(monkey).M1, ...
        ShenoyMonkeyData(monkey).Kinematics, ...
        kernelType,kernelScale,nonlinearL2,linearL2,holdOutNum, ...
        inferKinematicsFromLoadingsTrialAveraged(monkey).M1.nonlinearPrediction);
    % PMd.
    output(monkey).PMd = prepForAnalysis(ShenoyMonkeyData(monkey).PMd, ...
        ShenoyMonkeyData(monkey).Kinematics, ...
        kernelType,kernelScale,nonlinearL2,linearL2,holdOutNum, ...
        inferKinematicsFromLoadingsTrialAveraged(monkey).PMd.nonlinearPrediction);
end

end

%% Subfunction for organizing datasets before analysis.

function output = prepForAnalysis(data,kin,kernelType,kernelScale,nonLinearL2, ...
    linearL2,holdOutNum,nonlinearPreds)

% Prune out repeats.
data = pruneRepeats(data);
kin = pruneRepeats(kin);

% Assign the linear prediction.
output.nonlinearPrediction = nonlinearPreds;

%% Comparison against linear controls. 

% Perform linear decoding of position with 20 ms smoothing.
[lagged,output.linearDecoding.pos.kernel20ms.predicted] = predictKinUsingData(kin,holdOutNum,data, ...
    linearL2,10,0);
output.linearDecoding.pos.kernel20ms.varExplained = getVarExplainedKinematics( ...
    lagged,output.linearDecoding.pos.kernel20ms.predicted);
output.linearDecoding.pos.kernel20ms.posPVal = signrank( ...
    output.linearDecoding.pos.kernel20ms.varExplained.pos.array, ...
    output.nonlinearPrediction.varExplained.pos.array);
output.linearDecoding.pos.kernel20ms.velPVal = signrank( ...
    output.linearDecoding.pos.kernel20ms.varExplained.vel.array, ...
    output.nonlinearPrediction.varExplained.vel.array);

% Perform linear decoding of position with 100 ms smoothing.
[lagged,output.linearDecoding.pos.kernel100ms.predicted] = predictKinUsingData(kin,holdOutNum,data, ...
    linearL2,10,10);
output.linearDecoding.pos.kernel100ms.varExplained = getVarExplainedKinematics( ...
    lagged,output.linearDecoding.pos.kernel100ms.predicted);
output.linearDecoding.pos.kernel100ms.posPVal = signrank( ...
    output.linearDecoding.pos.kernel100ms.varExplained.pos.array, ...
    output.nonlinearPrediction.varExplained.pos.array);
output.linearDecoding.pos.kernel100ms.velPVal = signrank( ...
    output.linearDecoding.pos.kernel100ms.varExplained.vel.array, ...
    output.nonlinearPrediction.varExplained.vel.array);

% Perform linear decoding of velocity with 20 ms smoothing.
[lagged,output.linearDecoding.vel.kernel20ms.predicted] = predictKinUsingData(NthDerivative(kin,1,'raw'), ...
    holdOutNum,data,linearL2,10,0);
[output.linearDecoding.vel.kernel20ms.predicted,lagged] ...
    = postHoc(output.linearDecoding.vel.kernel20ms.predicted,lagged);
output.linearDecoding.vel.kernel20ms.varExplained = getVarExplainedKinematics( ...
    lagged,output.linearDecoding.vel.kernel20ms.predicted);
output.linearDecoding.vel.kernel20ms.posPVal = signrank( ...
    output.linearDecoding.vel.kernel20ms.varExplained.pos.array, ...
    output.nonlinearPrediction.varExplained.pos.array);
output.linearDecoding.vel.kernel20ms.velPVal = signrank( ...
    output.linearDecoding.vel.kernel20ms.varExplained.vel.array, ...
    output.nonlinearPrediction.varExplained.vel.array);

% Perform linear decoding of velocity with 100 ms smoothing.
[lagged,output.linearDecoding.vel.kernel100ms.predicted] = predictKinUsingData(NthDerivative(kin,1,'raw'), ...
    holdOutNum,data,linearL2,10,10);
[output.linearDecoding.vel.kernel100ms.predicted,lagged] = ...
    postHoc(output.linearDecoding.vel.kernel100ms.predicted,lagged);
output.linearDecoding.vel.kernel100ms.varExplained = getVarExplainedKinematics( ...
    lagged,output.linearDecoding.vel.kernel100ms.predicted);
output.linearDecoding.vel.kernel100ms.posPVal = signrank( ...
    output.linearDecoding.vel.kernel100ms.varExplained.pos.array, ...
    output.nonlinearPrediction.varExplained.pos.array);
output.linearDecoding.vel.kernel100ms.velPVal = signrank( ...
    output.linearDecoding.vel.kernel100ms.varExplained.vel.array, ...
    output.nonlinearPrediction.varExplained.vel.array);

%% Comparison against non-linear controls.

% Perform nonlinear decoding of position with 20 ms smoothing.
[lagged,output.nonlinearDecoding.pos.kernel20ms.predicted] = predictKinUsingDataNonlinear(kin,holdOutNum,data, ...
    kernelType,[kernelScale nonLinearL2],10,0);
output.nonlinearDecoding.pos.kernel20ms.varExplained = getVarExplainedKinematics( ...
    lagged,output.nonlinearDecoding.pos.kernel20ms.predicted);
output.nonlinearDecoding.pos.kernel20ms.posPVal = signrank( ...
    output.nonlinearDecoding.pos.kernel20ms.varExplained.pos.array, ...
    output.nonlinearPrediction.varExplained.pos.array);
output.nonlinearDecoding.pos.kernel20ms.velPVal = signrank( ...
    output.nonlinearDecoding.pos.kernel20ms.varExplained.vel.array, ...
    output.nonlinearPrediction.varExplained.vel.array);

% Perform nonlinear decoding of position with 100 ms smoothing.
[lagged,output.nonlinearDecoding.pos.kernel100ms.predicted] = predictKinUsingDataNonlinear(kin,holdOutNum,data, ...
    kernelType,[kernelScale nonLinearL2],10,10);
output.nonlinearDecoding.pos.kernel100ms.varExplained = getVarExplainedKinematics( ...
    lagged,output.nonlinearDecoding.pos.kernel100ms.predicted);
output.nonlinearDecoding.pos.kernel100ms.posPVal = signrank( ...
    output.nonlinearDecoding.pos.kernel100ms.varExplained.pos.array, ...
    output.nonlinearPrediction.varExplained.pos.array);
output.nonlinearDecoding.pos.kernel100ms.velPVal = signrank( ...
    output.nonlinearDecoding.pos.kernel100ms.varExplained.vel.array, ...
    output.nonlinearPrediction.varExplained.vel.array);

% Perform nonlinear decoding of velocity with 20 ms smoothing.
[lagged,output.nonlinearDecoding.vel.kernel20ms.predicted] = ...
    predictKinUsingDataNonlinear(NthDerivative(kin,1,'raw'),holdOutNum,data, ...
    kernelType,[kernelScale nonLinearL2],10,0);
[output.nonlinearDecoding.vel.kernel20ms.predicted,lagged] ...
    = postHoc(output.nonlinearDecoding.vel.kernel20ms.predicted,lagged);
output.nonlinearDecoding.vel.kernel20ms.varExplained = getVarExplainedKinematics( ...
    lagged,output.nonlinearDecoding.vel.kernel20ms.predicted);
output.nonlinearDecoding.vel.kernel20ms.posPVal = signrank( ...
    output.nonlinearDecoding.vel.kernel20ms.varExplained.pos.array, ...
    output.nonlinearPrediction.varExplained.pos.array);
output.nonlinearDecoding.vel.kernel20ms.velPVal = signrank( ...
    output.nonlinearDecoding.vel.kernel20ms.varExplained.vel.array, ...
    output.nonlinearPrediction.varExplained.vel.array);

% Perform nonlinear decoding of velocity with 100 ms smoothing.
[lagged,output.nonlinearDecoding.vel.kernel100ms.predicted] = ...
    predictKinUsingDataNonlinear(NthDerivative(kin,1,'raw'),holdOutNum,data, ...
    kernelType,[kernelScale nonLinearL2],10,10);
[output.nonlinearDecoding.vel.kernel100ms.predicted,lagged] ...
    = postHoc(output.nonlinearDecoding.vel.kernel100ms.predicted,lagged);
output.nonlinearDecoding.vel.kernel100ms.varExplained = getVarExplainedKinematics( ...
    lagged,output.nonlinearDecoding.vel.kernel100ms.predicted);
output.nonlinearDecoding.vel.kernel100ms.posPVal = signrank( ...
    output.nonlinearDecoding.vel.kernel100ms.varExplained.pos.array, ...
    output.nonlinearPrediction.varExplained.pos.array);
output.nonlinearDecoding.vel.kernel100ms.velPVal = signrank( ...
    output.nonlinearDecoding.vel.kernel100ms.varExplained.vel.array, ...
    output.nonlinearPrediction.varExplained.vel.array);

end



