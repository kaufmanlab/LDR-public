function [kinematics,predictedKin] = predictKinUsingData(kinematics,leaveOut,data,L2,lag,smoothing)
%% OVERVIEW

% This function attempts to predict held-out kinematics from their
% corresponding recordings using linear regression to fit a map from
% kinematics instantaneously to neural activity using linear regression. 

% 'The parameter L2 specifies an L2 penality on the regressors. 

% Smoothing specifies the Gaussian kernel.

%% Predict the held-out data.

% Prepare the data by lagging it and preparing the kinematic feature
% predictors. 
for cond = 1:size(data,2)
    if smoothing > 2
        % Get the new std for a kernel that updates the kernel to the
        % correct width. 
        newStd = (smoothing^2-2^2)^0.5;
        for neuron = 1:size(data(cond).matrix,1)
        data(cond).matrix(neuron,:) = ...
            imgaussfilt(data(cond).matrix(neuron,:),newStd, ...
            'padding','Symmetric','filterSize',101);
        end
    end
    data(cond).matrix = [data(cond).matrix(:,1:end-lag); ...
        ones(1,size(data(cond).matrix(:,1:end-lag),2))];
end
predictedKin = kinematics;
for cond = 1:size(data,2)
    kinematics(cond).matrix = [kinematics(cond).X kinematics(cond).Y];
    kinematics(cond).matrix = kinematics(cond).matrix(1+lag:end,:).';
end

% Create partitions of the data.
[train,test] = splitTrials([data().condNum],leaveOut);

% Loop through the partition.
for partition = 1:size(train,1)
    % Each col is vector in R^NC.
    heldInData = [data(train(partition,:)).matrix];
    % Get the held-in kinematics.
    heldInKin = [kinematics(train(partition,:)).matrix];
    linearMap = heldInKin*heldInData.'*pinv(heldInData*heldInData.' ...
        +L2*eye(size(heldInData,1)));
    % For the held out conditions, predict activity.
    predConds = find(test(partition,:));
    for cond = 1:length(predConds)
        predictedKin(predConds(cond)).matrix = ...
            linearMap*data(predConds(cond)).matrix;
    end
end
for cond = 1:size(data,2)
    predictedKin(cond).matrix = predictedKin(cond).matrix.';
    kinematics(cond).matrix = kinematics(cond).matrix.';
end
for cond = 1:size(data,2)
    predictedKin(cond).position = predictedKin(cond).matrix;
     predictedKin(cond).velocity = [0 0; diff(predictedKin(cond).position)];
end

end