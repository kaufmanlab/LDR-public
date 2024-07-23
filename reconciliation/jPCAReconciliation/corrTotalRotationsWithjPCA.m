function [corrCoeffs,pVals,jPCAVar,loadingsVar] ...
    = corrTotalRotationsWithjPCA(loadings,basisFxns,params,rotations)
%% OVERVIEW

% This function examines the correlation between total rotational activity,
% as expressed by the oscillatory temporal basis functions, and the
% magnitude of rotations in the jPCA plane.

%% Examine correlations.

% Assign two arrays for the covariates.
loadingsVar = zeros(size(loadings,2),1);
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
    loadingsVar(cond) = sum(nanvar(rotations(cond).softNorms.'.* ...
        ((loadings(cond).matrix(:,inds)-meanL(:,inds)) ...
        *basisFxns(temporalIndices,inds).'),0,2));
    for plane = 1:size(rotations(1).matrix,1)/2
        jPCAVar(cond,plane) = sum(nanvar( ...
            rotations(cond).matrix([1:2]+(plane-1)*2,:),0,2));
    end
end

% Assign the outputs.
[corrCoeffs,pVals] = corr(jPCAVar,loadingsVar);

end