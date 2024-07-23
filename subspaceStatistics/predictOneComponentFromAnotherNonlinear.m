function [groundTruthTensor,predictedData,reassambleTensor] = predictOneComponentFromAnotherNonlinear( ...
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
    condInds = find(train(partition,:));
    % Get loading matrices.
    [basisFxns,loadings,params] = eigTransform(data(train(partition,:)),factorDim);
    % Loop over the indices to get subcomponents of the manifold. 
    indices = [1 1+kron(1:floor(factorDim/2),[1 1])];
    for predictorIndex = 1:ceil(factorDim/2)
        for predictedIndex = 1:ceil(factorDim/2)
            % Get the inidices to predict from.
            predictor = (indices == predictorIndex);
            predicted = (indices == predictedIndex);
            % Assemble the covariates.
            predictedLoadings = loadings;
            predictorLoadings = loadings;
            for cond = 1:size(loadings,2)
                predictedLoadings(cond).matrix = ...
                    predictedLoadings(cond).matrix(:,predicted);
                predictedLoadings(cond).matrix = ...
                    predictedLoadings(cond).matrix(:);
                groundTruthTensor(predictedIndex,condInds(cond)).matrix ...
                    = loadings(cond).matrix(:,predicted)*basisFxns(:,predicted).';
                reassambleTensor(predictedIndex,condInds(cond)).matrix ...
                    = loadings(cond).matrix(:,predicted)*params.rawBasisFxns(:,predicted).';
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
                predictedData(predictorIndex,predictedIndex,predConds(cond)).matrix = ...
                    reshape(interpMap*interpCoeffs.', ...
                    size(data(1).matrix,1),length(find(predicted)))*basisFxns(:,predicted).';
            end
        end
    end
end

end