function trials = generateTrials(conds)
%% OVERVIEW

% This functions takes in conditions' activities in bins of
% ten-milliseconds, interpolates the firing rates to single millisecond
% precision, and then generates trials as a Poisson process. 

%% Parameters.

% Number of trials per condition.
trialCount = 24;

%% Generate trials.

% For each condition, generate trials.
for cond = 1:size(conds,2)
    activity = conds(cond).matrix;
    activity = interp1(-350:10:850,activity.',-350:850).';
    for trial = 1:trialCount
        trials(trial+(cond-1)*trialCount).matrix = poissrnd(activity/1000);
        trials(trial+(cond-1)*trialCount).condNum = conds(cond).condNum;
    end
end


end