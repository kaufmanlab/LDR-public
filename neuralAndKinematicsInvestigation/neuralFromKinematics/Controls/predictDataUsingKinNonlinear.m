function [data,predictedData] = predictDataUsingKinNonlinear(data,leaveOut,kinematics,LV,kernelType,params,lag,holdOutType)
%% OVERVIEW

% This function attempts to predict held-out neural recordings from their
% corresponding kinematics using kernel regression. There are a few options
% provided for how to parameterize the kinematic state at any given time. 

% 'The params vectors specifies 1) the length scale of the kernel and 2)
% the magnitude of the ridge regression.

%% Predict the held-out data.

% Prepare the data by lagging it and preparing the kinematic feature
% predictors. 

% Prepare a full model parameterization.
if strcmp(LV,'Pos')
    for cond = 1:size(data,2)
        data(cond).matrix = data(cond).matrix(:,1:end-lag);
    end
    predictedData = data;
    for cond = 1:size(data,2)
        kinematics(cond).matrix = ...
            [kinematics(cond).X kinematics(cond).Y ...
            (kinematics(cond).Y.^2 + kinematics(cond).X.^2).^0.5];
        kinematics(cond).matrix = kinematics(cond).matrix(1+lag:end,:).';
    end
elseif strcmp(LV,'Vel')
    for cond = 1:size(data,2)
        data(cond).matrix = data(cond).matrix(:,2:end-lag);
    end
    predictedData = data;
    for cond = 1:size(data,2)
        kinematics(cond).matrix = ...
            [imgaussfilt(diff(kinematics(cond).X),2, ...
            'Padding','symmetric','filterSize',21) imgaussfilt(diff(kinematics(cond).Y),2, ...
            'Padding','symmetric','filterSize',21) ...
            (imgaussfilt(diff(kinematics(cond).Y),2, ...
            'Padding','symmetric','filterSize',21).^2 + imgaussfilt(diff(kinematics(cond).X),2, ...
            'Padding','symmetric','filterSize',21).^2).^0.5];
        kinematics(cond).matrix = kinematics(cond).matrix(1+lag:end,:).';
    end
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
    heldInData = heldInData(:,1:2:end);
    % Get the held-in kinematics.
    % Get the kernel.
    heldInKin = [kinematics(train(partition,:)).matrix];
    heldInKin = heldInKin(:,1:2:end);
    kernelMat = evalKernel(heldInKin,heldInKin,kernelType,params(1));
    % Get the interpolation map.
    interpMap = heldInData*pinv(kernelMat+params(2)*eye(size(heldInKin,2)));
    % For the held out conditions, predict activity.
    predConds = find(test(partition,:));
    for cond = 1:length(predConds)
        % Get the interpolation coefficients.
        interpCoeffs = evalKernel(kinematics(predConds(cond)).matrix, ...
            heldInKin,kernelType,params(1));
        predictedData(predConds(cond)).matrix = ...
            interpMap*interpCoeffs.';
    end
end

end