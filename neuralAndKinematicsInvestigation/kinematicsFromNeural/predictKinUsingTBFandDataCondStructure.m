function [kinematics,predictedKinematics] = ...
    predictKinUsingTBFandDataCondStructure(kinematics,leaveOut,data,L2,factorDim)
%% OVERVIEW

% This function attempts to predict kinematics from the corresponding
% neural activity, in particular using the loading matrix of TBF as 
% regressors. This attempts to linearly link the dynamics of motor cortex
% with produced reaches. This function uses linear regression to link the
% two. 

% The parameter L2 specifies an L2 penality on the regressors. 

%% Predict the held-out data.

% Prepare the kinematics.
predictedKinematics = kinematics;
kinFxns = kinFuncEval(-350:10:850,'position');
for cond = 1:size(kinematics,2)
    kinematics(cond).matrix = [kinematics(cond).X kinematics(cond).Y];
    kinematics(cond).matrix = ...
        pinv(kinFxns)*kinematics(cond).matrix;
    kinematics(cond).matrix = kinematics(cond).matrix(:);
end

% Create partitions of the data.
[train,test] = splitTrials([data().condNum],leaveOut);

% Loop through the partition.
for partition = 1:size(train,1)
    % Get loading matrices.
    [basisFxns,loadings] = getBasisFxns(data(train(partition,:)),factorDim);
    for cond = 1:size(loadings,2)
        loadings(cond).matrix = loadings(cond).matrix(:);
    end
    heldInData = [loadings().matrix];
    heldInData = [heldInData; ones(1,size(heldInData,2))];
    % Get the held-in kinematics.
    heldInKin = [kinematics(train(partition,:)).matrix];
    % Get the map.
    linearMap = heldInKin*heldInData.'*pinv(heldInData*heldInData.' ...
        +L2*eye(size(heldInData,1)));
    % Test new conditions.
    predConds = find(test(partition,:));
    for cond = 1:length(predConds)
        % Get the loading matrix for the held-out condition.
        testLoading = data(predConds(cond)).matrix*pinv(basisFxns.');
        testLoading = testLoading(:);
        % Get the prediction of the data.
        predictedKinematics(predConds(cond)).matrix = ...
            linearMap(:,1:end-1)*testLoading+linearMap(:,end);
    end
end

% Correct for compression of the kinematics.
kinFxnsVel = kinFuncEval(-350:10:850,'velocity');
for cond = 1:size(kinematics,2)
    predictedKinematics(cond).velocity = kinFxnsVel*reshape( ...
        predictedKinematics(cond).matrix, ...
        [2 2]);
    predictedKinematics(cond).position = kinFxns*reshape( ...
        predictedKinematics(cond).matrix, ...
        [2 2]);
    kinematics(cond).matrix = [kinematics(cond).X kinematics(cond).Y];
end

end