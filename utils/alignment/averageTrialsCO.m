function trialAv = averageTrialsCO(trials,smoothKernel)
%% OVERVIEW

% This is a function for parcellating reaches to targets into components.

%% Parcellate trials.

% Classify trials by condition.
condInfo = [vertcat(trials.kinematics.targetTo) vertcat(trials.kinematics.targetFrom)];
inOut = [trials.kinematics.out];
condNum = 1;
useInds = 1:size(condInfo,1);
while ~isempty(condInfo)
    dists = squareform(pdist(condInfo));
    dists = dists(:,1);
    inds = find(dists == 0);
    for ind = inds.'
        trials.neuralActivity(useInds(ind)).condNum = condNum;
    end
    conds(condNum).out = inOut(useInds(inds(1)));
    if inOut(useInds(inds(1)))
        conds(condNum).target = trials.kinematics(useInds(inds(1))).targetTo;
    else
        conds(condNum).target = trials.kinematics(useInds(inds(1))).targetFrom;
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

% Loop over conditions, find "pair" condition.
for cond = find([conds.out] == 1)
    useTarg = conds(cond).target;
    inds = find([conds.out] == 0);
    matchTarg = vertcat(conds(inds).target);
    dists = sum((useTarg-matchTarg).^2,2);
    conds(cond).matchCond = inds(find(dists == 0));
end

% Loop over conditions.
indexes = find([conds.out] == 1);
for cond = 1:length(find([conds.out] == 1))
    % Get the trial indices for the outward section.
    trialInds = find([trials.neuralActivity.condNum] == indexes(cond));
    % Get the minimum end of the reaches, for the outward section.
    reachEnd = min([trials.kinematics(trialInds).end]);
    extentUse = min(extent(trialInds));
    outMove = unique([trials.kinematics(trialInds).moveInds]);
    % Assign initial arrays for neural activity.
    neuralOut = 0;
    kinOut = 0;
    for trial = trialInds
        neuralOut = neuralOut + 50*trials.neuralActivity(trial).matrix(:,1:reachEnd+extentUse)/length(trialInds);
        kinOut = kinOut + trials.kinematics(trial).matrix(1:reachEnd+extentUse,:)/length(trialInds);
    end
    % Get the trial indices for the inward section.
    trialInds = find([trials.neuralActivity.condNum] == conds(indexes(cond)).matchCond);
    % Get the minimum end of the reaches, for the inward section.
    reachEnd = min([trials.kinematics(trialInds).end]);
    extentUse = min(extent(trialInds));
    inMove = size(neuralOut,2)+unique([trials.kinematics(trialInds).moveInds]);
    % Assign initial arrays for neural activity.
    neuralIn = 0;
    kinIn = 0;
    for trial = trialInds
        neuralIn = neuralIn + 50*trials.neuralActivity(trial).matrix(:,1:reachEnd+extentUse)/length(trialInds);
        kinIn = kinIn + trials.kinematics(trial).matrix(1:reachEnd+extentUse,:)/length(trialInds);
    end
    % Concatenate.
    neural = [neuralOut neuralIn];
    kin = [kinOut; kinIn];
    if smoothKernel > 0
        % Smooth the neural activity.
        for channel = 1:size(neural,1)
            neural(channel,:) = imgaussfilt(neural(channel,:),smoothKernel, ...
                'padding','symmetric','filterSize',ceil((smoothKernel*10)/2)*2+1);
        end
    end
    % Assign everything.
    trialAv.neuralActivity(cond).matrix = neural;
    trialAv.kinematics(cond).matrix= kin;
    trialAv.kinematics(cond).targetOut = trials.kinematics(trialInds(1)).targetFrom;
    trialAv.kinematics(cond).moveInds = ismember(1:size(neural,2),unique([outMove inMove]));
    trialAv.kinematics(cond).angle = atan2(trialAv.kinematics(cond).targetOut(3), ...
        trialAv.kinematics(cond).targetOut(2));
end

end