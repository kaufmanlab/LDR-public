function recoverability = getRecoverability(data,loadings,basisFxns)
%% OVERVIEW

% This function takes data, loadings, and basis functions and returns the
% approximations of the basis functions recovered from the data, along with
% the variance explained.

%% Quantify the variance explained.

% Remove parts of the data not approximated, assuming that approximation is
% uniform across neurons but not time.
inds = find(isnan(basisFxns(:,1)));
basisFxns(inds,:) = [];
for cond = 1:size(data,2)
    data(cond).matrix(:,inds) = [];
end

% Quantify the variance explained.
recoverability.array = zeros(1,size(data,2));
for cond = 1:size(data,2)
    recoverability.array(cond) = 1 -sum(var(pinv(loadings(cond).matrix)*data(cond).matrix - basisFxns.',0,2))/ ...
        sum(var(basisFxns.',0,2));
    recoverability.recovered(cond).matrix = pinv(loadings(cond).matrix)*data(cond).matrix;
end
recoverability.array(recoverability.array < 0) = 0;
recoverability.mean = mean(recoverability.array);
recoverability.std = std(recoverability.array);

end

