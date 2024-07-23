function [data,predictedData,condList] = predictDataUsingTBFandKinCondStructureNonlinear(data,leaveOut,kinematics, ...
    kernelType,params,factorDim,holdOutType)
%% OVERVIEW

% This function attempts to predict held-out neural recordings from their
% corresponding kinematics by exploting hypothesized shared structure
% between the loading matrices and kinematics. This means that essentially
% the dynamics of motor cortex are determined by the produced reach. This
% function therefore uses nonlinear regression to predict the loading
% matrix of a condition from that condition's kinematics in the form of
% kernel regression.

% 'The params vectors specifies 1) the length scale of the kernel and 2)
% the magnitude of the ridge regression.

% The kernelType specifies the type of kernel. 

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

condList = [];

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
    % Form the kernel matrix.
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
            reshape(interpMap*interpCoeffs.', ...
            reshapeParams)*basisFxns.';
        condList = [condList predConds(cond)];
    end
end

end























