function approx = idDynamicsSameEigVals(data,eigVals,eigVecs)
%% OVERVIEW

% This function takes in data locked to movement onset, the "average"
% eigenvalues found by analysis, and the condition-specific eigenvectors
% for each condition, and approximates each condition using variable
% eigenvectors, invariant eigenvalues.

%% Find the dynamics describing each conditions.

% For each condition, find dynamics. 
approx = data;
for cond = 1:size(data,2)
    % Assign dynamics to the data in the low-dimensional space. 
    dynMatrix = real(eigVecs(cond).eigVecs*diag(eigVals)*pinv(eigVecs(cond).eigVecs));
    % Assign the approximation.
    for t = 21:121
        approx(cond).matrix(:,t) = ...
            dynMatrix*approx(cond).matrix(:,t-1);
    end
    approx(cond).matrix(:,1:19) = approx(cond).matrix(:,1:19)*nan;
end
    
end

