function [basisFxns,loadings,params] = eigTransform(data,basisFxnNum)
%% OVERVIEW

% This function takes the raw basis functions identified by an SVD on the
% dataset and identifies the nearest basis functions described by linear
% dynamics, which have the advantage of having an analytic form. The
% loadings corresponding to the raw basis functions are re-oriented to
% account for the eigen-mode property of the new basis functions.

%% Find the closest dynamical modes to the basis functions.

% Find the raw basis functions.
[basisFxnsRaw,loadings] = getBasisFxns(data,basisFxnNum);

% Fit dynamics to the basis functions. 
dynMatrix = basisFxnsRaw(21:end,:).'*pinv(basisFxnsRaw(20:end-1,:).');

% Eigendecompose the dynamics matrix.
[~,eigVals] = eigs(dynMatrix,size(dynMatrix,1));
eigVals = bioEigs(diag(eigVals));

% Re-sort based on frequency.
[~,I] = sort(atan2(abs(imag(eigVals)),real(eigVals))/(2*pi));
eigVals = eigVals(I);

% Use dynamics to simulate basis functions.
basisFxnsSim = basisFxnsRaw;
for t = 21:(size(data(1).matrix,2)-20)
    basisFxnsSim(t,:) = (dynMatrix*basisFxnsSim(t-1,:).').';
end

% Identify the pairs in the basis functions.
ID = pairID(eigVals);

% Assign the new basis functions' parameters.
params.param = zeros(size(dynMatrix,1),3);
params.dyn = dynMatrix;
params.functionalForm = zeros(size(dynMatrix,1),3);
params.meaning = ...
    ["Frequency (cycles/sec)" "Phase (degrees)" "Half-life (sec)"];
for el = 1:length(unique(ID))
    inds = find(ID == el);
    useEig = eigVals(inds(1));
    [f,t,frequency,halflife] = extractParams(useEig);
    params.functionalForm(inds(1),:) = [f 0 t];
    params.param(inds(1),:) = [frequency 0 halflife];
    if length(inds) == 2
        params.functionalForm(inds(2),:) = [f -pi/2 t];
        params.param(inds(2),:) = [frequency -90 halflife];
    end
end

% Assign the new basis functions.
timeStamp = [nan(19,1); (0:(size(data(1).matrix,2)-20)).'];
basisFxnsParam = basisFxnsRaw;
for col = 1:size(basisFxnsRaw,2)
    basisFxnsParam(:,col) = ...
        exp(-params.functionalForm(col,3)*timeStamp) ...
        .*cos(2*pi*params.functionalForm(col,1)*timeStamp ...
        + params.functionalForm(col,2));
end

% Fit the loading matrices.
correction = pinv(basisFxnsSim(20:end,:))*basisFxnsParam(20:end,:);
for cond = 1:size(data,2)
    loadings(cond).matrix = ...
        loadings(cond).matrix*pinv(correction.');
    loadings(cond).correction = correction;
end

% Assign the basis.
basisFxns = basisFxnsParam;
params.rawBasisFxns = (correction.'*basisFxnsRaw.').';
    
end

