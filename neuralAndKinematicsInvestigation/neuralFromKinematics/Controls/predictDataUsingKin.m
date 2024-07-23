function [data,predictedData] = predictDataUsingKin(data,leaveOut,kinematics,L2,lag,holdOutType)
%% OVERVIEW

% This function attempts to predict held-out neural recordings from their
% corresponding kinematics using linear regression to fit a map from
% kinematics instantaneously to neural activity using linear regression. 

% 'The parameter L2 specifies an L2 penality on the regressors. 

%% Predict the held-out data.

% Prepare the data by lagging it and preparing the kinematic feature
% predictors. 
for cond = 1:size(data,2)
    data(cond).matrix = data(cond).matrix(:,3:end-lag);
end
predictedData = data;
for cond = 1:size(data,2)
    kinematics(cond).matrix = ...
        [ones(size(kinematics(cond).X,1)-2,1) ....
        kinematics(cond).X(3:end) kinematics(cond).Y(3:end) ...
        (kinematics(cond).Y(3:end).^2 + kinematics(cond).X(3:end).^2).^0.5 ...
        cos(atan2(kinematics(cond).Y(3:end),kinematics(cond).X(3:end))) ...
        sin(atan2(kinematics(cond).Y(3:end),kinematics(cond).X(3:end))) ...
        diff(kinematics(cond).X(2:end)) diff(kinematics(cond).Y(2:end)) ...
        (diff(kinematics(cond).Y(2:end)).^2 + diff(kinematics(cond).X(2:end)).^2).^0.5 ...
        imgaussfilt(diff(kinematics(cond).X,2),2) imgaussfilt(diff(kinematics(cond).Y,2),2) ...
        (imgaussfilt(diff(kinematics(cond).X,2),2).^2 + imgaussfilt(diff(kinematics(cond).Y,2),2).^2).^0.5];
    kinematics(cond).matrix = kinematics(cond).matrix(1+lag:end,:).';
end

% Create partitions of the data.
if strcmp(holdOutType,'rand')
    [train,test] = splitTrials([data().condNum],leaveOut);
elseif strcmp(holdOutType,'angles')
    [train,test] = splitTrialsAngles([data().condNum],leaveOut,kinematics);
end

% Loop through the partition.
for partition = 1:size(train,1)
    % Each col is vector in R^NC.
    heldInData = [data(train(partition,:)).matrix];
    % Get the held-in kinematics.
    heldInKin = [kinematics(train(partition,:)).matrix];
    linearMap = heldInData*heldInKin.'*pinv(heldInKin*heldInKin.' ...
        +L2*eye(size(heldInKin,1)));
    % For the held out conditions, predict activity.
    predConds = find(test(partition,:));
    for cond = 1:length(predConds)
        predictedData(predConds(cond)).matrix = ...
            linearMap*kinematics(predConds(cond)).matrix;
    end
end

end