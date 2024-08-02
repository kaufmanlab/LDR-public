function activity = assembleCond(trials)
%% OVERVIEW

% This function takes in single-trials and averages
% them and returns the trial-averaged activity.

%% Average the activity.

% Average the trials.
activity = 0;
for trial = 1:size(trials,2)
    activity = activity + trials(trial).matrix/size(trials,2);
end

end

