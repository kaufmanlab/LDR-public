function [train,test] = leaveTrialsOutTarget(trials,holdOut,folds)
%% OVERVIEW

% This function is used to hold out a specified number of trials of a
% certain type. 

%% Hold trials out. 

% Classify trials by condition.
condInfo = [vertcat(trials.targetTo) vertcat(trials.targetFrom)];
condNum = 1;
useInds = 1:size(condInfo,1);
condNums = zeros(1,size(condInfo,2));
while ~isempty(condInfo)
    dists = squareform(pdist(condInfo));
    dists = dists(:,1);
    inds = find(dists == 0);
    for ind = inds.'
        condNums(useInds(ind)) = condNum;
    end
    condNum = condNum + 1;
    condInfo(inds,:) = [];
    useInds(inds) = [];
end
condNum = condNum - 1;

% Assign train, test matrices.
train = zeros(folds,size(trials,2));
test = train;

% Loop over folds, hold one angle out. 
for fold = 1:folds
    % Loop over conditions.
    for cond = 1:condNum
        % Find the associated indices.
        inds = find(condNums == cond);
        % Scramble the indices.
        inds = inds(randperm(length(inds)));
        % Assign the indices.
        train(fold,inds(1:end-holdOut)) = 1;
        test(fold,inds(end-holdOut+1:end)) = 1;
    end
end
train = logical(train);
test = logical(test);

end