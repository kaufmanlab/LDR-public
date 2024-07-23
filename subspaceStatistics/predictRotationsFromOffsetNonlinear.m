function [data,predictedData] = predictRotationsFromOffsetNonlinear( ...
    data,leaveOut,L2,factorDim,svEnergy,kernelType,lengthScale)
%% OVERVIEW

% This function attempts to show that the entire activity manifold of motor
% cortex changes jointly by predicting one component of the activity
% manifold from another. 

%% Infer activity.

% Create partitions of the data.
[train,test] = splitTrials([data().condNum],leaveOut);

% Loop through the partition.
for partition = 1:size(train,1)
    % Get loading matrices.
    [basisFxns,loadings,params] = eigTransform(data(train(partition,:)),factorDim);
    % Get the inidices to predict from.
    predictor = 1;
    predicted = 2:size(basisFxns,2);
    % Assemble the covariates.
    predictedLoadings = loadings;
    predictorLoadings = loadings;
    for cond = 1:size(loadings,2)
        predictedLoadings(cond).matrix = ...
            predictedLoadings(cond).matrix(:,predicted);
        predictedLoadings(cond).matrix = ...
            predictedLoadings(cond).matrix(:);
        predictorLoadings(cond).matrix = ...
            predictorLoadings(cond).matrix(:,predictor);
        predictorLoadings(cond).matrix = ...
            predictorLoadings(cond).matrix(:);
    end
    heldInPredictor = [predictorLoadings().matrix];
    heldInPredicted = [predictedLoadings().matrix];
    % Reduce the dimensionablity of the representation.
    [sv,S,~] = svd(heldInPredictor,'econ');
    S = diag(S);
    S = cumsum(S)/sum(S);
    indKeep = find(S > svEnergy);
    sv = sv(:,1:indKeep(1));
    heldInPredictor = sv.'*heldInPredictor;
    % Form the kernel matrix.
    kernelMat = evalKernel(heldInPredictor,heldInPredictor,kernelType,lengthScale);
    % Get the interpolation map.
    interpMap = heldInPredicted*pinv(kernelMat+L2*eye(size(heldInPredictor,2)));
    % Test new conditions.
    predConds = find(test(partition,:));
    for cond = 1:length(predConds)
        % Get the loading matrix for the held-out condition.
        testLoading = data(predConds(cond)).matrix*pinv(params.rawBasisFxns.');
        testLoading = testLoading(:,predictor);
        testCond = sv.'*testLoading(:);
        % Get the interpolation coefficients.
        interpCoeffs = evalKernel(testCond, ...
            heldInPredictor,kernelType,lengthScale);
        % Get the prediction of the data.
        predictedData(predConds(cond)).matrix = ...
            [testLoading(:,predictor) ...
            reshape(interpMap*interpCoeffs.', ...
            size(data(1).matrix,1),length(find(predicted)))]*basisFxns.';
    end
end

end