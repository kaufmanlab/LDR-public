function nullData = produceShuffledFPDataset(data,factorDim)
%% OVERVIEW

% This function is for controlling the mapping from the tuning vector/fixed
% point to the rotations. To control for this, a null dataset is needed.
% This is done by TBF-ing a dataset, then for each condition taking the
% projection of that dataset as a "proxy" for the temporal basis functions,
% then choosing two new conditions randomly, one of which controls the
% rotaitonal subspace and the other of which defines the tuning vector. 

%% Assign the dataset.

% Perform factorization.
[~,loadings,~] = eigTransform(data,factorDim);

% Shuffle the dataset.
nullData = data;
for cond = 1:size(nullData,2)
    condTuning = randi([1 size(nullData,2)],1,1);
    condRotations = randi([1 size(nullData,2)],1,1);
    subspace = [loadings(condTuning).matrix(:,1) loadings(condRotations).matrix(:,2:end)];
    nullData(cond).matrix = subspace*pinv(loadings(cond).matrix)*nullData(cond).matrix;
end

end