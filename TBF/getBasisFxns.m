function [basisFxns,loadings] = getBasisFxns(data,basisFxnNum)
%% OVERVIEW

% This function uses an SVD to extract the raw basis functions from a
% dataset. This is done by concatenating the datasets into an NC-by-T
% matrix, performing an SVD, then using the top k temporal singular vectors
% as the raw basis function. To avoid overfitting, the basis functions are
% extracted from a dataset where firing rates are normalized.

%% Get basis functions.

% Normalize firing rates.
dataUse = data;
maxFRs = max([data().matrix],[],2);
for cond = 1:size(data,2)
    dataUse(cond).matrix = dataUse(cond).matrix./maxFRs;
end

% Get the basis functions. 
[~,~,basisFxns] = svd(vertcat(data().matrix),'econ');
basisFxns = basisFxns(:,1:basisFxnNum);

% Get the loadings.
loadings = data;
for cond = 1:size(data,2)
    loadings(cond).matrix = data(cond).matrix*pinv(basisFxns.');
end

end

