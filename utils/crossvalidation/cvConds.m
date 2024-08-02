function [trainConds,testConds] = cvConds(trials,ratio)
%% OVERVIEW

% This function takes in trials and returns independent estimations of a
% condition's firing rates by averaging independent trials
% from the same condition.

%% Split and average trials.

% Get the conditions to assemble.
conds = unique([trials().condNum]);

% For each condition, split the corresponding trials (approximately) in
% half then pre-process and average.
for cond = conds
    inds = find([trials().condNum] == cond);
    inds = inds(randperm(length(inds)));
    trainInds = inds(1:round(end*ratio));
    testInds = inds(1+round(end*ratio):end);
    condInd = find(conds == cond);
    trainConds(condInd).matrix = assembleCond(trials(trainInds));
    testConds(condInd).matrix = assembleCond(trials(testInds));
end

end

