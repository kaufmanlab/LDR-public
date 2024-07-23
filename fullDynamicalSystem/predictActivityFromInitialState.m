function [data,predictedData] = predictActivityFromInitialState( ...
    data,leaveOut,initialStates,L2,factorDim,svEnergy)
%% OVERVIEW

% This function attempts to show that the initial state of neural activity
% is predictive of subsequent population states, arguing for a dynamical
% view of motor cortex activity. The initial state is quantified as the
% averaged firing rates 350-250 ms prior to movement onset, or 100 ms prior
% to movement modulation. The inference is bottlenecked through the loading
% matrices, which when expanded using the temporal basis functions provide
% an inference of neural activity. 

%% Infer activity.

% Create partitions of the data.
[train,test] = splitTrials([data().condNum],leaveOut);

% Loop through the partition.
for partition = 1:size(train,1)
    % Get loading matrices.
    [basisFxns,loadings,~] = eigTransform(data(train(partition,:)),factorDim);
    for cond = 1:size(loadings,2)
        loadings(cond).matrix = loadings(cond).matrix(:);
    end
    heldInData = [loadings().matrix];
    % Get the held-in kinematics.
    heldinInitialStates = initialStates(:,train(partition,:));
    [sv,S,~] = svd(heldinInitialStates,'econ');
    S = diag(S);
    S = cumsum(S)/sum(S);
    indKeep = find(S > svEnergy);
    sv = sv(:,1:indKeep(1));
    heldinInitialStates = [sv.'*heldinInitialStates; ...
        ones(1,size(heldinInitialStates,2))];
    % Get the map.
    linearMap = heldInData*heldinInitialStates.' ...
        *pinv(heldinInitialStates*heldinInitialStates.' ...
        +L2*eye(size(heldinInitialStates,1)));
    % Test new conditions.
    predConds = find(test(partition,:));
    for cond = 1:length(predConds)
        % Get the loading matrix for the held-out condition.
        testCond = sv.'*initialStates(:,predConds(cond));
        % Get the prediction of the data.
        predictedData(predConds(cond)).matrix = ...
            reshape(linearMap(:,1:end-1)*testCond+linearMap(:,end), ...
            size(data(1).matrix,1),factorDim)*basisFxns.';
    end
end

end