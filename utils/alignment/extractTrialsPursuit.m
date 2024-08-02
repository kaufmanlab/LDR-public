function [trials,trialAv] = extractTrialsPursuit(data,timePeriod,smoothKernel)
%% OVERVIEW

% This is a function for parcellating reaches to targets into components.

%% Parcellate trials.

% Correct kinematic offset.
data.Kinematics.ActualPos(:,2) = data.Kinematics.ActualPos(:,2) + 0.06;
data.TaskStateMasks.target(2,:) = data.TaskStateMasks.target(2,:) + 0.06;

% Get blocks where presentation is occuring. 
blocksOut = findBlocks(strcmp([data.TaskStateMasks.state_name],'Presentation'));

% Loop through, extract chunks specified.
count = 1;
for trial = 2:size(blocksOut,2)-1
    if blocksOut(trial+1).inds(1)+timePeriod(2)+1 <= size(data.SpikeCount,1) ...
        && blocksOut(trial).inds(end)-timePeriod(1)+1 > 0 ...
        && (blocksOut(trial).inds(1) - blocksOut(trial-1).inds(end)) < 500
        trials.neuralActivity(count).matrix = data.SpikeCount(...
            blocksOut(trial).inds(end)-timePeriod(1)+1:blocksOut(trial+1).inds(1)+timePeriod(2)+1, ...
            data.ActiveChannelMask).';
        trials.kinematics(count).matrix = data.Kinematics.ActualPos(...
            blocksOut(trial).inds(end)-timePeriod(1)+1:blocksOut(trial+1).inds(1)+timePeriod(2)+1,1:3);
        trials.targetOverTime(count).matrix = ...
            data.TaskStateMasks.target(1:3, ...
            blocksOut(trial).inds(end)-timePeriod(1)+1:blocksOut(trial+1).inds(1)+timePeriod(2)+1).';
        trials.kinematics(count).targetTo = ...
            data.TaskStateMasks.target(1:3,blocksOut(trial).inds(end)-5).';
        trials.kinematics(count).targetFrom = ...
            data.TaskStateMasks.target(1:3,blocksOut(trial).inds(1)-10).';
        trials.kinematics(count).end = ...
            blocksOut(trial+1).inds(1) + 10 - (blocksOut(trial).inds(end)-timePeriod(1)+1);
        plot(trials.kinematics(count).matrix(:,2),trials.kinematics(count).matrix(:,3))
        hold on
        scatter(trials.kinematics(count).targetTo(2),trials.kinematics(count).targetTo(3),50,[0 0 1])
        scatter(trials.kinematics(count).targetFrom(2),trials.kinematics(count).targetFrom(3),200,[1 0 0])
        xlim([-0.4 0.4])
        ylim([-0.4 0.4])
        pause(1/3)
        hold off
        count = count+1;
    end
end

% Smooth kinematics a bit.
for trial = 1:size(trials.kinematics,2)
    for dim = 1:3
        % Filter out nans.
        allInds = 1:size(trials.kinematics(trial).matrix,1);
        nanInds = isnan(trials.kinematics(trial).matrix(:,dim));
        goodInds = allInds;
        goodInds(nanInds) = [];
        trials.kinematics(trial).matrix(:,dim) = interp1(goodInds, ...
            trials.kinematics(trial).matrix(goodInds,dim),allInds, ...
            'nearest','extrap');
        % Smooth.
        trials.kinematics(trial).matrix(:,dim) = ...
            imgaussfilt(trials.kinematics(trial).matrix(:,dim),2,'padding', ...
            'symmetric','filterSize',31);
    end
end

% Classify trials by condition.
condInfo = [vertcat(trials.kinematics.targetTo) vertcat(trials.kinematics.targetFrom)];
condNum = 1;
useInds = 1:size(condInfo,1);
while ~isempty(condInfo) && ~isempty(pdist(condInfo))
    dists = squareform(pdist(condInfo));
    dists = dists(:,1);
    inds = find(dists == 0);
    for ind = inds.'
        trials.neuralActivity(useInds(ind)).condNum = condNum;
    end
    condNum = condNum + 1;
    condInfo(inds,:) = [];
    useInds(inds) = [];
end
condNum = condNum - 1;

% Find the minimum amount that trials extend past the end.
extent = zeros(size(trials.neuralActivity,2),1);
for trial = 1:size(trials.neuralActivity,2)
    extent(trial) = size(trials.neuralActivity(trial).matrix,2) ...
        - trials.kinematics(trial).end;
end
extent = floor(min(extent)/5)*5;

% Loop over conditions.
for cond = 1:condNum
    % Get the trial indices for the outward section.
    trialInds = find([trials.neuralActivity.condNum] == cond);
    % Get the minimum end of the reaches, for the outward section.
    reachEnd = min([trials.kinematics(trialInds).end]);
    % Assign initial arrays for neural activity.
    neural = 0;
    kin = 0;
    for trial = trialInds
        neural = neural + 50*trials.neuralActivity(trial).matrix(:,1:reachEnd+extent)/length(trialInds);
        kin = kin + trials.kinematics(trial).matrix(1:reachEnd+extent,:)/length(trialInds);
    end
    % Smooth the neural activity.
    for channel = 1:size(neural,1)
        neural(channel,:) = imgaussfilt(neural(channel,:),smoothKernel, ...
            'padding','symmetric','filterSize',ceil((smoothKernel*10)/2)*2+1);
    end
    % Assign everything.
    trialAv.neuralActivity(cond).matrix = neural;
    trialAv.kinematics(cond).matrix= kin;
    trialAv.kinematics(cond).targetOut = trials.kinematics(trialInds(1)).targetFrom;
end

end