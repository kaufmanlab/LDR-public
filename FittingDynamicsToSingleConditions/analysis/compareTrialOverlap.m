function output = compareTrialOverlap
%% OVERVIEW

% This function checks whether the variability in the subpace across
% condition is due to just estimation noise. To do this, the subspaces
% found by random partitions of trials are compared. 

%% Compare trials to the collected distribution.

% Load the data.
load('ShenoyMonkeyData');
ShenoyMonkeyData = ShenoyMonkeyData(1:2);
load('ShenoyMonkeyDataSingleTrial');

% Load the results.
load('condSpecificDynamics');

% For each monkey analyze.
for monkey = 1:size(ShenoyMonkeyData,2)
    output(monkey).M1 = analysis(condSpecificDynamics(monkey).M1.alignmentIndex, ....
        ShenoyMonkeyDataSingleTrial(monkey).M1, ...
        length(condSpecificDynamics(monkey).M1.averageEigVals), ...
        ShenoyMonkeyData(monkey).M1);
    output(monkey).PMd = analysis(condSpecificDynamics(monkey).PMd.alignmentIndex, ....
        ShenoyMonkeyDataSingleTrial(monkey).PMd, ...
        length(condSpecificDynamics(monkey).PMd.averageEigVals), ...
        ShenoyMonkeyData(monkey).PMd);
end

end

%% FUNCTION FOR ANALYZING DATA.

function output = analysis(estimateDistribution,trials,dim,data)

% Prune and preprocess the trials.
trials = preprocessTrials(pruneRepeatsTrial(trials));
data = pruneRepeats(data);

% Assign the number of partitions to use.
partitionNum = 500;

%% ALIGNMENT INDEX

% Assign the data and null distribution.
output.alignmentIndex.dataDist = estimateDistribution;
output.alignmentIndex.nullDist = zeros(length(unique([trials().condNum])),partitionNum);

% Estimate the null.
for partition = 1:partitionNum
    [conds1,conds2] = cvConds(trials,0.5);
    for cond = 1:size(conds1,2)
        [space,~,~] = svd(conds2(cond).matrix,'econ');
        output.alignmentIndex.nullDist(cond,partition) = ...
            getAlignmentIndex(conds1(cond).matrix,space(:,1:dim));
    end
end

% Test the significance.
useArray = output.alignmentIndex.dataDist;
for cond = 1:size(useArray,1)
    useArray(cond,cond) = nan;
end
useArray = useArray(:);
useArray(isnan(useArray)) = [];
output.alignmentIndex.histDist = useArray;
% Assign statistics.
output.alignmentIndex.dataDistMean = mean(useArray(:));
output.alignmentIndex.dataDistSTD = std(useArray(:));
output.alignmentIndex.nullDistMean = mean(output.alignmentIndex.nullDist(:));
output.alignmentIndex.nullDistSTD = std(output.alignmentIndex.nullDist(:));
% Assign the pVal.
output.alignmentIndex.pVal = ranksum(useArray,output.alignmentIndex.nullDist(:));
% Calculate ROC-AUC.
[~,~,~,output.alignmentIndex.AUC] = ...
    perfcurve([useArray*0; output.alignmentIndex.nullDist(:)*0+1], ...
    [useArray; output.alignmentIndex.nullDist(:)],1);
output.alignmentIndex.AUC = 0.5+abs(output.alignmentIndex.AUC-0.5);

% Calculate ROC-AUC.
[output.alignmentIndex.XCoord,output.alignmentIndex.YCoord,~,output.alignmentIndex.AUCMATLAB] = ...
    perfcurve([useArray*0; output.alignmentIndex.nullDist(:)*0+1], ...
    [useArray; output.alignmentIndex.nullDist(:)],1);
output.alignmentIndex.AUCMATLAB = 0.5+abs(output.alignmentIndex.AUCMATLAB-0.5);
for coord = 1:length(output.alignmentIndex.YCoord)
    if output.alignmentIndex.YCoord(coord) < output.alignmentIndex.XCoord(coord)
        storeXCoord = output.alignmentIndex.XCoord(coord);
        storeYCoord = output.alignmentIndex.YCoord(coord);
        output.alignmentIndex.XCoord(coord) = storeYCoord;
        output.alignmentIndex.YCoord(coord) = storeXCoord;
    end
end
output.alignmentIndex.AUCModified = sum((output.alignmentIndex.XCoord(2:end) ...
    - output.alignmentIndex.XCoord(1:end-1)).*output.alignmentIndex.YCoord(2:end));

end