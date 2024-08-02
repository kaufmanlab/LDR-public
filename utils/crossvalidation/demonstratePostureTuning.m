function output = demonstratePostureTuning(ShenoyMonkeyDataSingleTrial)
%% OVERVIEW

% This function is to demonstrate that rotations are not a trivial finding
% of neural activity starting and ending in the same place, because posture
% tuning means that post-reach activity is in a substantially-different
% location than before. This is done by computing firing rates (without
% smoothing), then comparing the distribution of inter-state distances
% between pre-movement modulation, post-movement, and between the two. 

%% Compute inter-state distances.

% Perform analysis for each dataset.
for monkey = 1:size(ShenoyMonkeyDataSingleTrial,2)
    output(monkey).M1 = analyzeStates(ShenoyMonkeyDataSingleTrial(monkey).M1);
    output(monkey).PMd = analyzeStates(ShenoyMonkeyDataSingleTrial(monkey).PMd);
end

end

%% Subfunction for performing analysis.

function out = analyzeStates(data)
% Get the unique conditions.
data = pruneRepeatsTrial(data);
conds = unique([data.condNum]);
% Loop over conditions.
for cond = 1:length(conds)
    % Assign activity.
    activity = 0;
    % Find the corresponding trials. 
    trials = find([data.condNum] == conds(cond));
    % Loop over trials, trial-average.
    for trial = trials
        activity = activity + rebin(full(double(data(trial).matrix)),10)/length(trials);
    end
    % Assign the activity. 
    out.neural(cond).matrix = activity;
    % Get the pairwise distance matrix. 
    distMat = squareform(pdist(activity.'));
    % Get the preparatory distances.
    prepDist = distMat(1:20,1:20);
    prepDist = prepDist(:);
    prepDist(prepDist == 0) = [];
    prepDist2(cond).vec = prepDist;
    % Get the preparatory distances.
    postDist = distMat(end-19:end,end-19:end);
    postDist = postDist(:);
    postDist(postDist == 0) = [];
    postDist2(cond).vec = postDist;
    % Get the inter-distances.
    interDist = distMat(1:20,end-19:end);
    interDist = interDist(:);
    interDist2(cond).vec = interDist;
end
% Concatenate.
out.interDist = [interDist2.vec];
out.postDist = [postDist2.vec];
out.prepDist = [prepDist2.vec];
out.pVals = [ranksum(out.interDist(:),out.postDist(:)) ...
    ranksum(out.interDist(:),out.prepDist(:))];
end