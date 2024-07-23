function output = compareTrialSimilarity
%% OVERVIEW

% This function checks whether the variability in the eigenvalues across
% condition is due to just estimation noise. To do this, the eigenvalues
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
    output(monkey).M1 = analysis(condSpecificDynamics(monkey).M1.eigValDist, ...
        ShenoyMonkeyDataSingleTrial(monkey).M1, ...
        length(condSpecificDynamics(monkey).M1.averageEigVals), ...
        condSpecificDynamics(monkey).M1.eigVals);
    output(monkey).PMd = analysis(condSpecificDynamics(monkey).PMd.eigValDist, ...
        ShenoyMonkeyDataSingleTrial(monkey).PMd, ...
        length(condSpecificDynamics(monkey).PMd.averageEigVals), ...
        condSpecificDynamics(monkey).PMd.eigVals);
end

end

%% FUNCTION FOR ANALYZING DATA.

function output = analysis(estimateDistributionEigVals,trials,dim,eigVals)

% Assign the number of partitions to use.
partitionNum = 500;

% Prune and preprocess the trials.
trials = preprocessTrials(pruneRepeatsTrial(trials));

% Assign the data and null distribution.
output.eigVals.dataDist = estimateDistributionEigVals;
output.eigVals.nullDist = zeros(length(unique([trials().condNum])),partitionNum);

% Assign an array for comparison eigenvalues.
compVals = zeros(size(eigVals,1),size(eigVals,2),partitionNum);

% Estimate the null.
for partition = 1:partitionNum
    [conds1,conds2] = cvConds(trials,0.5);
    % Fit dynamics to each partitions.
    eigVals1 = idDynamicsTruncated(conds1,dim);
    compVals(:,:,partition) = eigVals1;
    eigVals2 = idDynamicsTruncated(conds2,dim);
    output.eigVals.nullDist(:,partition) = ...
        eigValDistance(eigVals1,eigVals2);
end

% Assign the covariances for each eigenvalue.
for eVal = 1:size(compVals,1)
    for cond = 1:size(eigVals,2)
        if eVal >1
            useInds = find(abs(imag(compVals(eVal,cond,:))) > 0);
            output.params.cov(eVal,cond).matrix = cov([ ...
                    real(squeeze(compVals(eVal,cond,useInds))) imag(squeeze(compVals(eVal,cond,useInds)))]);
            output.params.mean(eVal,cond).vec = mean([ ...
                real(squeeze(compVals(eVal,cond,useInds))) imag(squeeze(compVals(eVal,cond,useInds)))]);
        end
    end
end

%% EIGVALS
% Test the significance.
useArray = output.eigVals.dataDist;
for cond = 1:size(useArray,1)
    useArray(cond,cond) = nan;
end
useArray = useArray(:);
useArray(isnan(useArray)) = [];
output.eigVals.histDist = useArray;
% Assign statistics.
output.eigVals.dataDistMean = mean(useArray(:));
output.eigVals.dataDistSTD = std(useArray(:));
output.eigVals.nullDistMean = mean(output.eigVals.nullDist(:));
output.eigVals.nullDistSTD = std(output.eigVals.nullDist(:));
% Assign the pVal.
output.eigVals.pVal = ranksum(useArray,output.eigVals.nullDist(:));
% Calculate ROC-AUC.
[output.eigVals.XCoord,output.eigVals.YCoord,~,output.eigVals.AUCMATLAB] = ...
    perfcurve([useArray*0; output.eigVals.nullDist(:)*0+1], ...
    [useArray; output.eigVals.nullDist(:)],1);
output.eigVals.AUCMATLAB = 0.5+abs(output.eigVals.AUCMATLAB-0.5);
for coord = 1:length(output.eigVals.YCoord)
    if output.eigVals.YCoord(coord) < output.eigVals.XCoord(coord)
        storeXCoord = output.eigVals.XCoord(coord);
        storeYCoord = output.eigVals.YCoord(coord);
        output.eigVals.XCoord(coord) = storeYCoord;
        output.eigVals.YCoord(coord) = storeXCoord;
    end
end
output.eigVals.AUCModified = sum((output.eigVals.XCoord(2:end) ...
    - output.eigVals.XCoord(1:end-1)).*output.eigVals.YCoord(2:end));

end

%% SUBFUNCTION FOR RAPIDLY FITTING DYNAMICS

function eigVals = idDynamicsTruncated(data,dimUse)

% For each condition, find dynamics. 
eigVals = zeros(dimUse,size(data,2));
for cond = 1:size(data,2)
    % Fit dynamics to the data. 
    dynMatrix = (data(cond).matrix(:,21:end)) ...
        *pinv(data(cond).matrix(:,20:end-1));
    % Solve RRR.
    [subspace,~,~] = svd(dynMatrix*data(cond).matrix(:,20:end),'econ');
    subspace = subspace(:,1:dimUse);
    dynMatrix = (subspace*subspace.')*dynMatrix;
    % Eigendecompose the dynamics.
    [~,e] = eigs(dynMatrix,dimUse);
    e = bioEigs(diag(e));
    [~,I] = sort(atan2(abs(imag(e)),real(e))/(2*pi));
    e = e(I);
    eigVals(:,cond) = e;
end
    
end

