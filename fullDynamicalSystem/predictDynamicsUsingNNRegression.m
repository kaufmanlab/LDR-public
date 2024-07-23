function [predictedActivity,data] = ...
    predictDynamicsUsingNNRegression(data,leaveOut,loadings,basisFxns)
%% OVERVIEW

% This function attempts to predict the evolution of motor cortex activity
% on a given condition from previously-observed conditions using NN
% regression. This is done by predicting the derivative of the state from
% the state itself. 

%% Predict the held-out data.

% Assign data.
temporalInds = find(~isnan(basisFxns(:,1)));
for cond = 1:size(data,2)
    data(cond).matrix = loadings(cond).matrix*basisFxns.';
end

% Declare the BETA parameter.
beta = 2;

% Prepare the output.
predictedActivity = data;

% Create partitions of the data.
[train,test] = splitTrials([data().condNum],leaveOut);

% Prepare the matrices to use for prediction and the state variables.
stateVars = data;
stateDiffs = data;
nanInds = find(isnan(basisFxns(:,1)));
for cond = 1:size(data,2)
    stateVars(cond).matrix = ...
        data(cond).matrix(:,temporalInds(1:end-1));
    stateDiffs(cond).matrix = ...
        diff(data(cond).matrix(:,temporalInds),1,2);
end

% Loop through the partition.
for partition = 1:size(train,1)
    % Get the held-in data.
    heldInData = [stateVars(train(partition,:)).matrix];
    heldInDiffs = [stateDiffs(train(partition,:)).matrix];
    % Test new conditions.
    predConds = find(test(partition,:));
    for cond = 1:length(predConds)
        % nan non-dynamical indices.
        predictedActivity(predConds(cond)).matrix(:,nanInds) = ...
            predictedActivity(predConds(cond)).matrix(:,nanInds)*nan;
        % Loop over timepoints.
        for timePoint = 21:121
            % Find the pairwise kernel similarity with the training
            % dataset.
            dists = pdist2( ...
                predictedActivity(predConds(cond)).matrix(:,timePoint-1).', ...
                heldInData.');
            dists = exp(-dists/beta)/sum(exp(-dists/beta));
            % Get the prediction of the derivative.
            predictedActivity(predConds(cond)).matrix(:,timePoint) = ...
                predictedActivity(predConds(cond)).matrix(:,timePoint-1) ...
                +heldInDiffs*dists.';
        end
    end
end

end