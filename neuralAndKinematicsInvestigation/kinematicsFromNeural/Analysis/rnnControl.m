function output = rnnControl
%% OVERVIEW

% This function compares the performance of linear and nonlinear prediction
% using TBF and prediction made a by standard RNN. 

%% Parameters.

% Declare the kernel used for non-linear mapping.
kernelType = 'OU';

% Declares the length scale of the kernel.
kernelScale = 10^3;

% Declare the regularization on the kernel for nonlinear mapping.
nonlinearL2 = 0.1;

% Declare the regularization on the linear mapping.
linearL2 = 1;

% Declare the dimensionality for decoding.
decodeDim = 20;

%% Compare the datasets for M1-N.

% Load the rnns inference. 
load('testRnnVarExplained');
output.rnnResults = testRnnVarExplained;

% Load the data.
load('ShenoyMonkeyDataSingleTrial');
ShenoyMonkeyDataSingleTrial(2).M1 = ...
    poolDatasets(ShenoyMonkeyDataSingleTrial(2).M1,ShenoyMonkeyDataSingleTrial(2).PMd);

% Assemble half the data into training and testing trials.
singleTrialUse = preprocessTrials(ShenoyMonkeyDataSingleTrial(2).M1(1:2:end));
singleTrialTest = ShenoyMonkeyDataSingleTrial(2).M1(2:2:end);

% Assemble test kinematics.
kinTest = ShenoyMonkeyDataSingleTrial(2).Kinematics(2:2:end);

% Assemble the trial averages.
kinTrain = ShenoyMonkeyDataSingleTrial(2).Kinematics(1:2:end);
uniqueConds = unique([singleTrialUse.condNum]);
for cond = 1:length(uniqueConds)
    useCond = uniqueConds(cond);
    inds = find([singleTrialUse.condNum] == useCond);
    useActivity(cond).matrix = assembleCond(singleTrialUse(inds));
    useActivity(cond).condNum = useCond;
    useKin(cond).X = 0;
    useKin(cond).Y = 0;
    for ind = inds
        useKin(cond).X = useKin(cond).X + kinTrain(ind).X/length(inds);
        useKin(cond).Y = useKin(cond).Y + kinTrain(ind).Y/length(inds);
    end
end

% Perform linear inference.
[output.data,output.linearPrediction.predicted] = ...
    predictKinUsingTBFandDataCondStructureSingletrial(useKin,kinTest,1, ...
    useActivity,singleTrialTest,linearL2,9);
output.linearPrediction.varExplained = getVarExplainedKinematics(output.data, ...
    output.linearPrediction.predicted);
output.linearPrediction.pVal = signrank(output.rnnResults, ...
    output.linearPrediction.varExplained.vel.array);

% Perform nonlinear inference.
[~,output.nonlinearPrediction.predicted] = ...
    predictKinUsingTBFandDataCondStructureNonlinearSingleTrial(useKin,kinTest,1,useActivity,singleTrialTest, ...
    kernelType,[kernelScale nonlinearL2 decodeDim],9);
output.nonlinearPrediction.varExplained = getVarExplainedKinematics(output.data, ...
    output.nonlinearPrediction.predicted);
output.nonlinearPrediction.pVal = signrank(output.rnnResults, ...
    output.nonlinearPrediction.varExplained.vel.array);

end