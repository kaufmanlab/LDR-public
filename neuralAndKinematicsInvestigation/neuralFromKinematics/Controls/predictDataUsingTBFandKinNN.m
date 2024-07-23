function [data,predictedData] = predictDataUsingTBFandKinNN(data,leaveOut,kinematics,factorDim,holdOutType)
%% OVERVIEW

% This function attempts to predict held-out neural recordings from their
% corresponding kinematics using nearest-neighbor regression. The variance
% explained by this model tests whether the kernel regression works
% simply by memorizing the training set. 

%% Predict the held-out data.

% Prepare the kinematics.
predictedData = data;
kinFxns = kinFuncEval(-350:10:850,'position');
for cond = 1:size(data,2)
    kinematics(cond).matrix = ...
        pinv(kinFxns)...
        *[kinematics(cond).X kinematics(cond).Y];
    kinematics(cond).matrix = kinematics(cond).matrix(:);
end

% Create partitions of the data.
if strcmp(holdOutType,'rand')
    [train,test] = splitTrials([data().condNum],leaveOut);
elseif strcmp(holdOutType,'angles')
    [train,test] = splitTrialsAngles([data().condNum],leaveOut,kinematics);
end

% Loop through the partition.
for partition = 1:size(train,1)
    % Get loading matrices.
    [basisFxns,loadings,~] = eigTransform(data(train(partition,:)),factorDim);
    reshapeParams = size(loadings(1).matrix);
    for cond = 1:size(loadings,2)
        loadings(cond).matrix = loadings(cond).matrix(:);
    end
    heldInData = [loadings().matrix];
    % Get the held-in kinematics.
    heldInKin = [kinematics(train(partition,:)).matrix];
    % For the held out conditions, predict activity.
    predConds = find(test(partition,:));
    for cond = 1:length(predConds)
        % Get the distances from the training set at the test point.
        dists = sum((kinematics(predConds(cond)).matrix - heldInKin).^2).^0.5;
        % Select the minimum.
        [~,I] = min(dists);
        % Return the nearest neighbor.
        predictedData(predConds(cond)).matrix = ...
            reshape(heldInData(:,I), ...
            reshapeParams)*basisFxns.';
    end
end

end























