function data = removePotent(data,EMG,kinematics)
%% OVERVIEW

% This function estimates the "potent" space of motor cortex by using
% linear regression to estimate the dimensions that encode for muscle
% activity, along with hand velocity and position. This is then rejected 
% from the data entirely. This is useful for showing that LDR-based 
% decoding does not rely on explicitly-encoded kinematic-related activity.

%% Estimate and reject potent space. 

% Copy data, truncate to account for differentiation of the hand position.
velData = data;
for cond = 1:size(data,2)
    velData(cond).matrix = velData(cond).matrix(:,1:end-1);
end

% Compile kinematics.
kinParams = zeros(size([velData.matrix],2),4);
for cond = 1:size(data,2)
    kinParams((1:size(velData(cond).matrix,2))+ ...
        (cond-1)*size(velData(cond).matrix,2),:) = ...
        [kinematics(cond).X(1:end-1) kinematics(cond).Y(1:end-1) ....
        diff(kinematics(cond).X) diff(kinematics(cond).Y)];
end

% Estimate the map from neural activity to kinematics.
kinParamsSub = (kinParams - mean(kinParams,2)).';
velDataSub = [velData.matrix]-mean([velData.matrix],2);
map = kinParamsSub*pinv(velDataSub);

% Loop over muscle groups.
for emgGroup = 1:size(EMG,2)
    % Get the associated conditions.
    useConds = EMG(emgGroup).numbers;
    % Compile neural data.
    neuralArray = [[data(useConds).prepare] [data(useConds).matrix]];
    neuralArray = neuralArray - mean(neuralArray,2);
    % Compile an array of muscle activity.
    muscleArray = [EMG(emgGroup).prepare.' EMG(emgGroup).activity.'];
    muscleArray = muscleArray - mean(muscleArray,2);
    % Append to the map.
    map = [map; muscleArray*pinv(neuralArray)];
end

% Loop over data, remove the projection on the encoding dimensions.
for cond = 1:size(data,2)
    data(cond).matrix = data(cond).matrix - pinv(map)*map*data(cond).matrix;
end

end