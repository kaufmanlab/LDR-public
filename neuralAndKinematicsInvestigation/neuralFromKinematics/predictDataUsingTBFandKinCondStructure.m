function [data,predictedData] = predictDataUsingTBFandKinCondStructure(data,leaveOut,kinematics,L2,factorDim,holdOutType)
%% OVERVIEW

% This function attempts to predict held-out neural recordings from their
% corresponding kinematics by exploting hypothesized shared structure
% between the loading matrices and kinematics. This means that essentially
% the dynamics of motor cortex are determined by the produced reach. THis
% function therefore uses linear regression to predict the loading
% matrix of a condition from that condition's kinematics.

% 'The parameter L2 specifies an L2 penality on the regressors. 

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
    heldInKin = [heldInKin; ones(1,size(heldInKin,2))];
    linearMap = heldInData*heldInKin.'*pinv(heldInKin*heldInKin.' ...
        +L2*eye(size(heldInKin,1)));
    % For the held out conditions, predict activity.
    predConds = find(test(partition,:));
    for cond = 1:length(predConds)
        predictedData(predConds(cond)).matrix = ...
            reshape(linearMap(:,1:end-1) ...
            *kinematics(predConds(cond)).matrix+linearMap(:,end), ...
            reshapeParams)*basisFxns.';
    end
end

end