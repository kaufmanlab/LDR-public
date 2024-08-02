function trialStruct = extractHumanTrials(data,keepObsOnly,customMask,pad,prePad)
%% OVERVIEW

% This function takes in a data structure produced by Sliman et al, in
% particular of a participant viewing a center-out reaching task. This is
% divided up into trials that are aligned to movement onset. 

%% Get trials.

% Find blocks that constitute trials. 
currentTrial = 0;
inReach = false;
inReach2 = false;
centerPresented = false;
sawEnd = false;
for timepoint = 1:length([data.TaskStateMasks.state_name])
    % Wait for presentation to begin tracking.
    if currentTrial == 0
        if strcmp(data.TaskStateMasks.state_name{timepoint},'Presentation1') 
            currentTrial = currentTrial+1;
            trial(currentTrial).inds = timepoint-prePad:timepoint;
        end
    % If a reach has already been seen, default to different instructions.
    else
        % If presentation and not currently in a reach, count it as part of
        % current reach.
        if strcmp(data.TaskStateMasks.state_name{timepoint},'Presentation1') && ~inReach
            trial(currentTrial).inds = [trial(currentTrial).inds timepoint];
        % If anything else is seen, count it as part of current reach. If a
        % reach is seen, count it as "in" a reach.
        elseif strcmp(data.TaskStateMasks.state_name{timepoint},'Reach') 
            trial(currentTrial).inds = [trial(currentTrial).inds timepoint];
            if ~inReach
                inReach = true;
                trial(currentTrial).moveStart = timepoint;
            end
        elseif strcmp(data.TaskStateMasks.state_name{timepoint},'Presentation2') 
            trial(currentTrial).inds = [trial(currentTrial).inds timepoint];
            if ~centerPresented
                trial(currentTrial).secondPresentation = timepoint;
                centerPresented = true;
            end
        elseif strcmp(data.TaskStateMasks.state_name{timepoint},'Center') 
            trial(currentTrial).inds = [trial(currentTrial).inds timepoint];
            if ~inReach2
                inReach2 = true;
                trial(currentTrial).move2Start = timepoint;
            end
            sawEnd = true;
        % If presentation is seen and current in a reach, then count up to
        % next trial and break out of reaching. Pad a bit to the end of
        % each trial.
        elseif (strcmp(data.TaskStateMasks.state_name{timepoint},'Presentation1') && inReach) ...
                || (strcmp(data.TaskStateMasks.state_name{timepoint},'InterTrial') && inReach)
            trial(currentTrial).inds = [trial(currentTrial).inds timepoint:timepoint+pad];
            if sawEnd
                trial(currentTrial).ended = true;
                sawEnd = false;
            else
                trial(currentTrial).ended = false;
            end
            currentTrial = currentTrial+1;
            trial(currentTrial).inds = timepoint-prePad:timepoint;
            % Reset every metric. 
            inReach = false;
            inReach2 = false;
            centerPresented = false;
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
            listStates(count).vec = ...
                data.TaskStateMasks.state_name(trial(thisTrial).inds);
            trialStruct.kinematics(count).targetOverTime = ...
                data.TaskStateMasks.target(1:3,trial(thisTrial).inds).';
            naninds = find(~isnan(sum(trialStruct.kinematics(count).matrix,2)));
            trialStruct.kinematics(count).matrix = interp1( ...
                naninds, ...
                trialStruct.kinematics(count).matrix(naninds,:), ...
                1:length(trialStruct.kinematics(count).matrix),'nearest','extrap');
            trialStruct.kinematics(count).target = data.TaskStateMasks.target(1:3,trial(thisTrial).moveStart);
            trialStruct.kinematics(count).angle = ...
                atan2(data.TaskStateMasks.target(2,trial(thisTrial).moveStart), ...
                data.TaskStateMasks.target(3,trial(thisTrial).moveStart));
            trialStruct.kinematics(count).move = trial(thisTrial).moveStart ...
                - trial(thisTrial).inds(1)+1;
            trialStruct.kinematics(count).move2 = trial(thisTrial).move2Start ...
                - trial(thisTrial).inds(1)+1;
            trialStruct.kinematics(count).secondPresentation = trial(thisTrial).secondPresentation ...
                - trial(thisTrial).inds(1)+1;
            trialStruct.kinematics(count).timeInds =  ...
                trial(thisTrial).inds-trial(thisTrial).inds(1);
            % Indicate if this was an observation trial or not.
            trialStruct.kinematics(count).observed = ...
                data.TaskStateMasks.brain_control_weight(1,trial(thisTrial).moveStart) == 0;
            trialStruct.kinematics(count).ended = ...
                trial(thisTrial).ended;
        end
    end
end

% If screen is selected, keep only observation trials.
allStates = unique(data.TaskStateMasks.state_name);
totalTrialNum = size(trialStruct.kinematics,2);
trialStruct.stateList = allStates;
if keepObsOnly
    %keepInds = intersect(find([trialStruct.kinematics.observed]),find([trialStruct.kinematics.ended]));
    keepInds = find([trialStruct.kinematics.ended]);
    trialStruct.kinematics = ...
        trialStruct.kinematics(keepInds);
    trialStruct.neuralActivity = ...
        trialStruct.neuralActivity(keepInds);
    % Time lock.
    minsOut = zeros(size(trialStruct.kinematics,2),1);
    maxesOut = zeros(size(trialStruct.kinematics,2),1);
    % Find the distribution of times, both for out and back.
    for trial = 1:size(trialStruct.kinematics,2)
        minsOut(trial) = ...
            trialStruct.kinematics(trial).move - 1;
        maxesOut(trial) = ...
            trialStruct.kinematics(trial).move2 - (trialStruct.kinematics(trial).move+1);
    end
    mostMinOut = min(minsOut);
    leastMaxOut = min(maxesOut);
    minsIn = zeros(size(trialStruct.kinematics,2),1);
    maxesIn = zeros(size(trialStruct.kinematics,2),1);
    % Find the distribution of times, both for out and back.
    for trial = 1:size(trialStruct.kinematics,2)
        minsIn(trial) = ...
            trialStruct.kinematics(trial).move2 - (trialStruct.kinematics(trial).move+leastMaxOut+1);
        maxesIn(trial) = ...
            trialStruct.kinematics(trial).timeInds(end) - (trialStruct.kinematics(trial).move2+1);
    end
    mostMinIn = min(minsIn);
    leastMaxIn = min(maxesIn);
    % Lock both the out and back portion of the trial.
    for trial = 1:size(trialStruct.kinematics,2)
        indsOut = (-mostMinOut:leastMaxOut)+trialStruct.kinematics(trial).move;
        indsIn = (-mostMinIn:leastMaxIn)+trialStruct.kinematics(trial).move2;
        allInds = union( ...
            indsOut,indsIn);
        useInds = find(ismember(trialStruct.kinematics(trial).timeInds,allInds));
        trialStruct.kinematics(trial).matrix = ...
            trialStruct.kinematics(trial).matrix(useInds,:);
        trialStruct.neuralActivity(trial).matrix = ...
            trialStruct.neuralActivity(trial).matrix(:,useInds);
        trialStruct.kinematics(trial).timeInds = ...
            trialStruct.kinematics(trial).timeInds(useInds);
        trialStruct.kinematics(trial).targetOverTime = ...
            trialStruct.kinematics(trial).targetOverTime(useInds,:);
        listStates(trial).vec = ...
            listStates(trial).vec(useInds);
        trialStruct.states(trial).list = zeros(length(allStates),length(listStates(trial).vec));
        for state = 1:length(allStates)
            trialStruct.states(trial).list(state,:) = ...
                ismember(listStates(trial).vec,allStates(state));
        end
        trialStruct.states(trial).stateIn = ...
            1+length(allStates)-sum(cumsum(trialStruct.states(trial).list));
    end
end

end