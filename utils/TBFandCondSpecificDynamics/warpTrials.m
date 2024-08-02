function [newConds,trialData,trialDuration] = warpTrials(trials,trialsExist,kinematicsCond)
%% OVERVIEW

% This function takes in trials, in the case of trialsExist = 1, or takes
% in data and generates trial otherwise, along with kinematic information,
% and returns trials that have been warped such that all reaches last the
% same duration. Trials are then pre-processed into new conditions. For
% plotting purposes the trial durations are also returned.

%% Warp trials into the same duration.

% Prepare trials.
if trialsExist == 1
    trialData = pruneRepeatsTrial(trials);
    %kinematics = pruneRepeatsTrial(kinematics);
else
    trialData = generateTrials(pruneRepeats(trials));
    %kinematicsCond = pruneRepeats(kinematics);
end
for trial = 1:size(trialData,2)
    kinematics(trial) = kinematicsCond(trialData(trial).condNum);
end

% Get the old and new duration of reaching, indexed to movement onset.
trialDuration = zeros(1,size(trialData,2));
for trial = 1:size(trialData,2)
    vel = (diff(kinematics(trial).X).^2+diff(kinematics(trial).Y).^2).^0.5;
    fastInds = find(vel >= max(vel)/10);
    % Quantify the end of the reach as the last moment the arm was moving
    % 10% of the maximum speed.
    trialDuration(trial) = 10*(fastInds(end)-36);
end

% For each trial, warp the trial to the average duration.
for trial = 1:size(trialData,2)
    trialData(trial).matrix = warpTrial(trialData(trial).matrix, ...
        trialDuration(trial),600);
end

% Pre-process the trial.
trialData = preprocessTrials(trialData);

% Assemble into conditions.
condInds = unique([trialData().condNum]);
for cond = 1:length(condInds)
    useTrials = [trialData().condNum] == condInds(cond);
    newConds(cond).matrix = assembleCond(trialData(useTrials));
    newConds(cond).condNum = condInds(cond);
    newConds(cond).oldDuration = mean(trialDuration(useTrials));
end

end