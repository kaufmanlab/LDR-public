function approx = useFactorization(loadings,basisFxns)
%% OVERVIEW

% This function takes loading matrices and basis functions and returns
% approximations to a dataset.

%% Produce the approximation.

% Create the approximation.
approx = loadings;
for cond = 1:size(loadings,2)
    approx(cond).matrix = approx(cond).matrix*basisFxns.';
end

end

