function output = kinToTBF
%% OVERVIEW

% As an analysis function, this function will load its own data. This
% function uses kinematics to guess the trial-averaged firing rates of the
% neural population, bottle-necked by guessing the loading matrix for the
% temporal basis functions of the condition. 

%% Parameters.

% Declare the kernel used for non-linear mapping.
kernelType = 'OU';

% Declares the length scale of the kernel.
kernelScale = 10^3;

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
load('basisFxnNum');

% For each monkey and brain region, prepare for analysis.
for monkey = 1:size(ShenoyMonkeyData,2)
    % M1.
    output(monkey).M1 = prepForAnalysis(ShenoyMonkeyData(monkey).M1, ...
        ShenoyMonkeyData(monkey).Kinematics, ...
        kernelType,kernelScale,nonlinearL2,linearL2,holdOutNum, ...
        basisFxnNum(monkey).M1.maxDim);
    % PMd.
    output(monkey).PMd = prepForAnalysis(ShenoyMonkeyData(monkey).PMd, ...
        ShenoyMonkeyData(monkey).Kinematics, ...
        kernelType,kernelScale,nonlinearL2,linearL2,holdOutNum, ...
        basisFxnNum(monkey).PMd.maxDim);
end

end

%% Subfunction for organizing datasets before analysis.

function output = prepForAnalysis(data,kin,kernelType,kernelScale,nonLinearL2, ...
    linearL2,holdOutNum,factorNum)

% Prune out repeats.
data = pruneRepeats(data);
kin = pruneRepeats(kin);

% Perform linear mapping with random hold-out.
[output.data,output.linearPrediction.predicted] = predictDataUsingTBFandKinCondStructure( ...
    data,holdOutNum,kin,linearL2,factorNum,'rand');
% Get the variance explained by this.
output.linearPrediction.varExplained = getVarExplained( ...
    output.data,output.linearPrediction.predicted,'ind');
% Get a sense of the variance explained after scrambling.
[~,null] = predictDataUsingTBFandKinCondStructure( ...
    data(randperm(size(data,2))),holdOutNum,kin,linearL2,factorNum,'rand');
output.linearPrediction.nullVarExplained = getVarExplained( ...
    output.data,null,'ind');
output.linearPrediction.pVal = signrank( ...
    output.linearPrediction.nullVarExplained.array, ...
    output.linearPrediction.varExplained.array);

% Perform nonlinear prediction with random hold-out.
[~,output.nonlinearPrediction.predicted] = predictDataUsingTBFandKinCondStructureNonlinear( ...
    data,holdOutNum,kin,kernelType,[kernelScale nonLinearL2],factorNum,'rand');
% Get the variance explained by this.
output.nonlinearPrediction.varExplained = getVarExplained( ...
    output.data,output.nonlinearPrediction.predicted,'ind');
% Get a sense of the variance explained after scrambling.
[~,null] = predictDataUsingTBFandKinCondStructureNonlinear( ...
    data(randperm(size(data,2))),holdOutNum,kin,kernelType, ...
    [kernelScale nonLinearL2],factorNum,'rand');
output.nonlinearPrediction.nullVarExplained = getVarExplained( ...
    output.data,null,'ind');
output.nonlinearPrediction.rand.pVal = signrank( ...
    output.nonlinearPrediction.nullVarExplained.array, ...
    output.nonlinearPrediction.varExplained.array);

% Perform nonlinear prediction with structured hold-out.
[~,output.nonlinearPrediction.control.predicted] = predictDataUsingTBFandKinCondStructureNonlinear( ...
    data,holdOutNum,kin,kernelType,[kernelScale nonLinearL2],factorNum,'angles');
% Get the variance explained by this.
output.nonlinearPrediction.control.varExplained = getVarExplained( ...
    output.data,output.nonlinearPrediction.control.predicted,'ind');
end

















