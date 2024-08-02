function trialStruct = extractHumanTrialsTargetPursuit(data,keepObsOnly,customMask)
%% OVERVIEW

% This function takes in a data structure produced by Sliman et al, in
% particular of a participant viewing a center-out reaching task. This is
% divided up into trials that are aligned to movement onset. 

%% Get trials.

% Find blocks that constitute trials. 
currentTrial = 0;
inReach = false;
sawEnd = false;
for timepoint = 1:length([data.TaskStateMasks.state_name])
    % Wait for presentation to begin tracking.
    if currentTrial == 0
        if strcmp(data.TaskStateMasks.state_name{timepoint},'Presentation') 
            currentTrial = currentTrial+1;
            trial(currentTrial).inds = timepoint;
        end
    % If a reach has already been seen, default to different instructions.
    else
        % If presentation and not currently in a reach, count it as part of
        % current reach.
        if strcmp(data.TaskStateMasks.state_name{timepoint},'Presentation') && ~inReach
            trial(currentTrial).inds = [trial(currentTrial).inds timepoint];
        % If anything else is seen, count it as part of current reach. If a
        % reach is seen, count it as "in" a reach.
        elseif strcmp(data.TaskStateMasks.state_name{timepoint},'Reach') 
            trial(currentTrial).inds = [trial(currentTrial).inds timepoint];
            if ~inReach
                inReach = true;
                trial(currentTrial).moveStart = timepoint;
            end
        % If presentation is seen and current in a reach, then count up to
        % next trial and break out of reaching. Pad a bit to the end of
        % each trial.
        elseif strcmp(data.TaskStateMasks.state_name{timepoint},'Presentation') && inReach
            trial(currentTrial).inds = [trial(currentTrial).inds timepoint:timepoint+15];
            currentTrial = currentTrial+1;
            trial(currentTrial).inds = timepoint;
            inReach = false;
        end
    end
end

% Loop over trials, grab trials.
count = 0;
for thisTrial = 1:size(trial,2)
    if size(data.SpikeCount,1)  > trial(thisTrial).inds(end)
        if ~isnan(prod(prod(data.SpikeCount(trial(thisTrial).inds,data.ActiveChannelMask).')))
            count = count+1;
            % Get neural activity.
            if isempty(customMask)
                trialStruct.neuralActivity(count).matrix = ...
                    data.SpikeCount(trial(thisTrial).inds,data.ActiveChannelMask).';
            else
                trialStruct.neuralActivity(count).matrix = ...
                    data.SpikeCount(trial(thisTrial).inds,customMask).';
            end
            % Grab kinematics, target.
            trialStruct.kinematics(count).matrix = ...
                data.Kinematics.ActualPos(trial(thisTrial).inds,1:3);
            trialStruct.states(count).list = ...
                data.TaskStateMasks.state_name(trial(thisTrial).inds);
            naninds = find(~isnan(sum(trialStruct.kinematics(count).matrix,2)));
            trialStruct.kinematics(count).matrix = interp1( ...
                naninds, ...
                trialStruct.kinematics(count).matrix(naninds,:), ...
                1:length(trialStruct.kinematics(count).matrix),'nearest','extrap');
            trialStruct.kinematics(count).target = data.TaskStateMasks.target(1:3,trial(thisTrial).moveStart);
            trialStruct.kinematics(count).from = trialStruct.kinematics(count).matrix(1,:).';
            trialStruct.kinematics(count).angle = ...
                atan2(data.TaskStateMasks.target(2,trial(thisTrial).moveStart), ...
                data.TaskStateMasks.target(3,trial(thisTrial).moveStart));
            trialStruct.kinematics(count).move = trial(thisTrial).moveStart ...
                - trial(thisTrial).inds(1)+1;
            trialStruct.kinematics(count).timeInds =  ...
                trial(thisTrial).inds - trial(thisTrial).moveStart;
            % Indicate if this was an observation trial or not.
            trialStruct.kinematics(count).observed = ...
                data.TaskStateMasks.brain_control_weight(1,trial(thisTrial).moveStart) == 0;
        end
    end
end

% If screen is selected, keep only observation trials.
if keepObsOnly
    keepInds = find([trialStruct.kinematics.observed]);
    %keepInds = find([trialStruct.kinematics.observed]);
    trialStruct.kinematics = ...
        trialStruct.kinematics(keepInds);
    trialStruct.neuralActivity = ...
        trialStruct.neuralActivity(keepInds);
    % Time lock.
    mins = zeros(size(trialStruct.kinematics,2),1);
    maxes = zeros(size(trialStruct.kinematics,2),1);
    for trial = 1:size(trialStruct.kinematics,2)
        mins(trial) = min(trialStruct.kinematics(trial).timeInds);
        maxes(trial) = max(trialStruct.kinematics(trial).timeInds);
    end
    mostMin = max(mins);
    leastMax = min(maxes);
    for trial = 1:size(trialStruct.kinematics,2)
        useInds = find(ismember(trialStruct.kinematics(trial).timeInds,mostMin:leastMax));
        trialStruct.kinematics(trial).matrix = ...
            trialStruct.kinematics(trial).matrix(useInds,:);
        trialStruct.neuralActivity(trial).matrix = ...
            trialStruct.neuralActivity(trial).matrix(:,useInds);
        trialStruct.kinematics(trial).timeInds = ...
            trialStruct.kinematics(trial).timeInds(useInds);
        trialStruct.states(trial).list = ...
            trialStruct.states(trial).list(useInds);
    end
end

end