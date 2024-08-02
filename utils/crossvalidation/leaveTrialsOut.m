function [train,test] = leaveTrialsOut(angles,holdOut,folds)
%% OVERVIEW

% This function is used to hold out a specified number of trials of a
% certain type. 

%% Hold trials out. 

% Get the unique type of each trial.
conds = unique(angles);

% Assign train, test matrices.
train = zeros(folds,length(angles));
test = train;

% Loop over folds, hold one angle out. 
for fold = 1:folds
    % Loop over conditions.
    for cond = conds
        % Find the associated indices.
        inds = find(angles == cond);
        % Scramble the indices.
        inds = inds(randperm(length(inds)));
        % Assign the indices.
        train(fold,inds(1:end-holdOut)) = 1;
        test(fold,inds(end-holdOut+1:end)) = 1;
    end
end
train = logical(train);
test = logical(test);

end