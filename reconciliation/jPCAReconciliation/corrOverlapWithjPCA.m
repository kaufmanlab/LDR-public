function [corrCoeffs,pVals,jPCAVar,overlap] ...
    = corrOverlapWithjPCA(loadings,basisFxns,params,rotations)
%% OVERVIEW

% This function examines the correlation between the overlap between
% rotational activity and the jPCA planes, and the variance within the jPCA
% planes.

%% Examine correlations.

% Assign two arrays for the covariates.
overlap = zeros(size(loadings,2),size(rotations(1).matrix,1)/2);
jPCAVar = zeros(size(loadings,2),size(rotations(1).matrix,1)/2);

% Get the rotational indices.
inds = find(params(:,1) > 0);
temporalIndices = find(~isnan(basisFxns(:,1)));

% Get the meanL.
meanL = 0;
for cond = 1:size(loadings,2)
    meanL = meanL + loadings(cond).matrix/size(loadings,2);
end

% For each condition, quantify covariates.
for cond = 1:size(loadings,2)
    for plane = 1:size(rotations(1).matrix,1)/2
        overlap(cond,plane) = getAlignmentIndex(...
            rotations(cond).softNorms.'.* ...
            ((loadings(cond).matrix(:,inds)-meanL(:,inds)) ...
            *basisFxns(temporalIndices,inds).'), ...
            rotations(cond).jPCs(:,[1:2]+(plane-1)*2));
        jPCAVar(cond,plane) = sum(nanvar( ...
            rotations(cond).matrix([1:2]+(plane-1)*2,:),0,2));
    end
end

% Assign the outputs.
[corrCoeffs,pVals] = corr(jPCAVar,overlap);
corrCoeffs = diag(corrCoeffs);
pVals = diag(pVals);

end